require 'strscan'
require 'thread'

module Bio
  module MAF

    class ParseError < Exception; end

    class Header
      attr_accessor :vars, :alignment_params

      def initialize(vars, params)
        @vars = vars
        @alignment_params = params
      end

      def version
        vars[:version]
      end

      def scoring
        vars[:scoring]
      end

    end

    class Block
      attr_reader :vars, :sequences, :offset, :size

      def initialize(*args)
        @vars, @sequences, @offset, @size = args
      end

      def raw_seq(i)
        sequences.fetch(i)
      end

      def each_raw_seq
        sequences.each { |s| yield s }
      end

      def text_size
        sequences.first.text.size
      end

    end

    class Sequence
      attr_reader :source, :start, :size, :strand, :src_size, :text
      attr_accessor :i_data, :quality
      alias_method :source_size, :src_size

      def initialize(*args)
        @source, @start, @size, @strand, @src_size, @text = args
      end

      def write_fasta(writer)
        writer.write("#{source}:#{start}-#{start + size}",
                     text)
      end
    end

    class ChunkReader
      attr_accessor :chunk_size, :chunk_shift, :pos
      attr_reader :f
      def initialize(f, chunk_size)
        @f = f
        self.chunk_size = chunk_size
        @pos = 0
      end

      def chunk_size=(size)
        check_chunk_size(size)
        @chunk_size = size
        # power of 2 so don't worry about rounding
        @chunk_shift = Math.log2(size).to_i
      end

      def check_chunk_size(size)
        if size < 1
          raise "Invalid chunk size: #{size}"
        end
        ## test whether it is a power of 2
        ## cf. http://bit.ly/JExNc4
        if size & (size - 1) != 0
          raise "Invalid chunk size (not a power of 2): #{size}}"
        end
      end

      # Reads the next chunk of the file.
      def read_chunk
        chunk = f.read(@chunk_size)
        @pos += chunk.bytesize if chunk
        return chunk
      end

     def read_chunk_at(offset, size_hint=@chunk_size)
        f.seek(offset)
        chunk = f.read(size_hint)
        @pos = offset + chunk.bytesize
        return chunk
      end
   end

   class ThreadedChunkReader < ChunkReader
     attr_reader :f

     def initialize(f, chunk_size, buffer_size=64)
       super(f, chunk_size)
       @buffer = SizedQueue.new(buffer_size)
       @eof_reached = false
       start_read_ahead
     end

     def start_read_ahead
       @read_thread = Thread.new { read_ahead }
     end

     def read_ahead
       # n = 0
       begin
         f_pos = 0
         until f.eof?
           chunk = f.read(@chunk_size)
           @buffer << [f_pos, chunk]
           f_pos += chunk.bytesize
           # n += 1
           # if (n % 100) == 0
           #   $stderr.puts "buffer size: #{@buffer.size}"
           # end
         end
         @eof_reached = true
       rescue Exception
         @read_ahead_ex = $!
         $stderr.puts "read_ahead aborting: #{$!}"
       end
     end

     def read_chunk
       raise "readahead failed: #{@read_ahead_ex}" if @read_ahead_ex
       if @eof_reached && @buffer.empty?
         return nil
       else
         c_pos, chunk = @buffer.shift()
         @pos = c_pos
         return chunk
       end
     end

   end

   class Parser

      ## Parses alignment blocks by reading a chunk of the file at a time.

      attr_reader :header, :file_spec, :f, :s, :cr, :at_end
      attr_reader :chunk_start, :last_block_pos
      attr_accessor :sequence_filter

      SEQ_CHUNK_SIZE = 131072
      RANDOM_CHUNK_SIZE = 4096
      MERGE_MAX = SEQ_CHUNK_SIZE

      def initialize(file_spec, opts={})
        chunk_size = opts[:chunk_size] || SEQ_CHUNK_SIZE
        @random_access_chunk_size = opts[:random_chunk_size] || RANDOM_CHUNK_SIZE
        @merge_max = opts[:merge_max] || MERGE_MAX
        @chunk_start = 0
        @file_spec = file_spec
        @f = File.open(file_spec)
        reader = opts[:chunk_reader] || ChunkReader
        @cr = reader.new(@f, chunk_size)
        @s = StringScanner.new(read_chunk())
        set_last_block_pos!
        @at_end = false
        _parse_header()
      end

      BLOCK_START = /^(?=a)/
      BLOCK_START_OR_EOS = /(?:^(?=a))|\z/
      EOL_OR_EOF = /\n|\z/

      def read_chunk
        cr.read_chunk
      end

      def fetch_blocks(fetch_list, filters=nil)
        ## fetch_list: array of [offset, length, block_count] tuples
        ## returns array of Blocks
        return fetch_blocks_merged(merge_fetch_list(fetch_list))
      end

      def fetch_blocks_merged(fetch_list)
        r = []
        old_chunk_size = cr.chunk_size
        cr.chunk_size = @random_access_chunk_size
        @at_end = false
        begin
          # inhibit fragment joining
          @last_block_pos = -1
          fetch_list.each do |offset, len, block_offsets|
            if (chunk_start <= offset) \
              && (offset < (chunk_start + s.string.size))
              ## the selected offset is in the current chunk
              s.pos = offset - chunk_start
            else
              chunk = cr.read_chunk_at(offset, len)
              @chunk_start = offset
              @s = StringScanner.new(chunk)
            end
            # read chunks until we have the entire merged set of
            # blocks ready to parse
            while s.string.size < len
              s.string << cr.read_chunk()
            end
            # parse the blocks
            block_offsets.each do |expected_offset|
              block = parse_block
              unless block
                raise "Parse error: expected a block at offset #{expected_offset} but could not parse one!"
              end
              unless block.offset == expected_offset
                raise "Parse error: got block with offset #{block.offset}, expected #{expected_offset}!"
              end
              r << block
            end
          end
          return r
        ensure
          cr.chunk_size = old_chunk_size
          set_last_block_pos!
        end
      end

      def merge_fetch_list(orig_fl)
        fl = orig_fl.dup
        r = []
        until fl.empty? do
          cur = fl.shift
          if r.last \
            && (r.last[0] + r.last[1]) == cur[0] \
            && (r.last[1] + cur[1]) <= @merge_max
            # contiguous with the previous one
            # add to length and increment count
            r.last[1] += cur[1]
            r.last[2] << cur[0]
          else
            cur << [cur[0]]
            r << cur
          end
        end
        return r
      end

      def set_last_block_pos!
        @last_block_pos = s.string.rindex(BLOCK_START)
      end

      def _parse_header
        parse_error("not a MAF file") unless s.scan(/##maf\s*/)
        vars = parse_maf_vars()
        align_params = nil
        while s.scan(/^#\s*(.+?)\n/)
          if align_params == nil
            align_params = s[1]
          else
            align_params << ' ' << s[1]
          end
        end
        @header = Header.new(vars, align_params)
        s.skip_until BLOCK_START || parse_error("Cannot find block start!")
      end

      ## On finding the start of a block:
      ## See whether we are at the last block in the chunk.
      ##   If at the last block:
      ##     If at EOF: last block.
      ##     If not:
      ##       Read the next chunk
      ##       Find the start of the next block in that chunk
      ##       Concatenate the two block fragments
      ##       Parse the resulting block
      ##       Promote the next scanner, positioned

      def trailing_nl(string)
        if string.empty?
          false
        else
          s.string[s.string.size - 1] == "\n"
        end
      end

      def parse_block
        return nil if at_end
        if s.pos != last_block_pos
          # in non-trailing block
          parse_block_data
        else
          # in trailing block fragment
          parse_trailing_fragment
        end
      end

      ## Read chunks and accumulate a leading fragment until we
      ## encounter a block start or EOF.
      def gather_leading_fragment
        leading_frag = ''
        while true
          next_chunk_start = cr.pos
          next_chunk = read_chunk
          if next_chunk
            next_scanner = StringScanner.new(next_chunk)
            # If this trailing fragment ends with a newline, then an
            # 'a' at the beginning of the leading fragment is the
            # start of the next alignment block.
            if trailing_nl(leading_frag) || trailing_nl(s.string)
              pat = BLOCK_START
            else
              pat = /(?:\n(?=a))/
            end
            frag = next_scanner.scan_until(pat)
            if frag
              # got block start
              leading_frag << frag
              break
            else
              # no block start in this
              leading_frag << next_chunk
            end
          else
            # EOF
            @at_end = true
            break
          end
        end
        return leading_frag, next_scanner, next_chunk_start
      end

      def parse_trailing_fragment
        leading_frag, next_scanner, next_chunk_start = gather_leading_fragment
        # join fragments and parse
        trailing_frag = s.rest
        joined_block = trailing_frag + leading_frag
        @chunk_start = chunk_start + s.pos
        @s = StringScanner.new(joined_block)
        begin
          block = parse_block_data
        rescue ParseError => pe
          parse_error "Could not parse joined fragments: #{pe}\nTRAILING: #{trailing_frag}\nLEADING: #{leading_frag}"
        end
        # Set up to parse the next block
        @s = next_scanner
        @chunk_start = next_chunk_start
        unless @at_end
          set_last_block_pos!
        end
        return block
      end

      def parse_error(msg)
        s_start = [s.pos - 10, 0].max
        s_end = [s.pos + 10, s.string.length].min
        if s_start > 0
          left = s.string[s_start..(s.pos - 1)]
        else
          left = ''
        end
        right = s.string[s.pos..s_end]
        extra = "pos #{s.pos} [#{chunk_start + s.pos}], last #{last_block_pos}"

        raise ParseError, "#{msg} at: '#{left}>><<#{right}' (#{extra})"
      end

      S = 's'.getbyte(0)
      I = 'i'.getbyte(0)
      E = 'e'.getbyte(0)
      Q = 'q'.getbyte(0)
      COMMENT = '#'.getbyte(0)

      def parse_block_data
        block_start_pos = s.pos
        block_offset = chunk_start + block_start_pos
        s.scan(/^a\s*/) || parse_error("bad a line")
        block_vars = parse_maf_vars()
        seqs = []
        payload = s.scan_until(/^(?=a)/)
        unless payload
          payload = s.rest
          s.pos = s.string.size # jump to EOS
        end
        lines = payload.split("\n")
        until lines.empty?
          line = lines.shift
          first = line.getbyte(0)
          if first == S
            seq = parse_seq_line(line)
            seqs << seq if seq
          # when 'i'
          #   parts = line.split
          #   parse_error("wrong i source #{src}!") unless seqs.last.source == src
          #   seqs.last.i_data = parts.slice(2..6)
          # when 'q'
          #   _, src, quality = line.split
          #   parse_error("wrong q source #{src}!") unless seqs.last.source == src
          #   seqs.last.quality = quality
          elsif [I, E, Q, COMMENT, nil].include? first
            next
          else
            parse_error "unexpected line: '#{line}'"
          end
        end
        return Block.new(block_vars,
                         seqs,
                         block_offset,
                         s.pos - block_start_pos)
      end

      def parse_seq_line(line)
        _, src, start, size, strand, src_size, text = line.split
        return nil if sequence_filter && ! seq_filter_ok?(src)
        begin
          Sequence.new(src,
                       start.to_i,
                       size.to_i,
                       STRAND_SYM.fetch(strand),
                       src_size.to_i,
                       text)
        rescue KeyError
          parse_error "invalid sequence line: #{line}"
        end
      end

      def seq_filter_ok?(src)
        if sequence_filter[:only_species]
          src_sp = src.split('.', 2)[0]
          m = sequence_filter[:only_species].find { |sp| src_sp == sp }
          return m
        else
          return true
        end
      end

      def parse_maf_vars
        vars = {}
        while s.scan(/(\w+)=(\S*)\s+/) do
          vars[s[1].to_sym] = s[2]
        end
        vars
      end

      STRAND_SYM = {
        '+' => :+,
        '-' => :-
      }

      def each_block
        until at_end
          yield parse_block()
        end
      end

    end

  end
  
end
