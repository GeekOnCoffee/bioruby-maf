require 'kyotocabinet'
require 'bitstring'
require 'bindata'

require 'bio-ucsc-api'
require 'bio-genomic-interval'

module Bio

  module MAF

    class IndexKey < BinData::Record
      endian :big

      uint8    :marker, :value => 0xFF
      uint8    :seq_id
      uint16   :bin
      uint32   :seq_start
      uint32   :seq_end
    end

    class IndexValue < BinData::Record
      endian :big

      uint64   :b_offset
      uint32   :b_length
      uint32   :text_size
      uint8    :n_seq
      uint64   :species_vec
    end

    module KVHelpers

      KEY_FMT = "CCS>L>L>"
      KEY_SCAN_FMT = "xCS>L>L>"
      CHROM_BIN_PREFIX_FMT = "CCS>"
      VAL_FMT = "Q>L>L>CQ>"
      VAL_IDX_OFFSET_FMT = "Q>L>"
      VAL_TEXT_SIZE_FMT = "@12L>"
      VAL_SPECIES_FMT = "@17Q>"

      module_function

      def extract_species_vec(entry)
        entry[1].species_vec
      end

      def extract_index_offset(entry)
        [entry[1].b_offset, entry[1].b_length]
      end

      def extract_text_size(entry)
        entry[1].text_size
      end

      def bin_start_prefix(seq_id, bin)
        k = IndexKey.new
        k.seq_id = seq_id
        k.bin = bin
        return k.to_binary_s.slice(0...(k.bin.offset + k.bin.num_bytes))
      end
    end

    class KyotoIndex
      include KVHelpers

      attr_reader :db, :species
      attr_accessor :index_sequences

      MAX_SPECIES = 64

      ## Key-value store index format
      ##
      ## This format is designed for Kyoto Cabinet but should work on
      ## other key-value databases allowing binary data.
      ##
      ## Index metadata is stored as ASCII text, but index data is
      ## stored as packed binary values.
      ##
      ## Index metadata:
      ##
      ##   Sequence IDs:
      ##     sequence:<name> => <id>
      ##
      ##     Each indexed sequence has a corresponding entry of this
      ##     kind. The <name> parameter is the sequence or chromosome
      ##     name as found in the MAF file, e.g. mm8.chr7. The <id>
      ##     parameter is assigned when the sequence is indexed, and
      ##     can be from 0 to 255.
      ##
      ## Index data:
      ##
      ##   For each sequence upon which an index is built, one index
      ##   entry is generated per MAF alignment block. The key
      ##   identifies the sequence, the UCSC index bin, and the
      ##   zero-based start and end positions of the sequence. The
      ##   value gives the offset and size of the alignment block
      ##   within the MAF file.
      ##
      ##   All values are stored as big-endian, unsigned packed binary
      ##   data.
      ##
      ## Keys: (12 bytes) [CCS>L>L>]
      ##
      ##   0xFF (1 byte):
      ##      index entry prefix
      ##   Sequence chromosome ID (1 byte):
      ##      corresponds to sequence:<name> entries
      ##   UCSC bin (16 bits)
      ##   Sequence start, zero-based, inclusive (32 bits)
      ##   Sequence end, zero-based, exclusive (32 bits)
      ##
      ## Values (25 bytes) [Q>L>L>CQ>]
      ##
      ##   MAF file offset (64 bits)
      ##   MAF alignment block length (32 bits)
      ##   Block text size (32 bits)
      ##   Number of sequences in block (8 bits)
      ##   Species bit vector (64 bits)
      ##
      ## Example:
      ##
      ##  For a block with sequence 0, bin 1195, start 80082334, end
      ##       80082368, MAF offset 16, and MAF block length 1087:
      ##
      ##     |  |id| bin | seq_start | seq_end   |
      ## key: FF 00 04 AB 04 C5 F5 9E 04 C5 F5 C0
      ##
      ##     |         offset        |  length   |   ts   |ns|  species_vec  |
      ## val: 00 00 00 00 00 00 00 10 00 00 04 3F  [TODO]

      #### Public API

      # Open an existing index for reading.
      def self.open(path)
        return KyotoIndex.new(path)
      end

      # Build a new index from the MAF file being parsed by PARSER,
      # and store it in PATH.
      def self.build(parser, path)
        idx = self.new(path)
        idx.build_default(parser)
        return idx
      end

      # Find all alignment blocks in the genomic regions in the list
      # of Bio::GenomicInterval objects INTERVALS, and parse them with
      # PARSER.
      def find(intervals, parser, filter={})
        parser.fetch_blocks(fetch_list(intervals, filter))
      end

      # Close the underlying Kyoto Cabinet database handle.
      def close
        db.close
      end

      #### KyotoIndex Internals

      def initialize(path)
        @key = IndexKey.new
        @val = IndexValue.new
        @species = {}
        @species_max_id = -1
        if (path.size > 1) and File.exist?(path)
          mode = KyotoCabinet::DB::OREADER
        else
          mode = KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE
        end
        @db = KyotoCabinet::DB.new
        @path = path
        unless db.open(path, mode)
          raise "Could not open DB file!"
        end
        if mode == KyotoCabinet::DB::OREADER
          load_index_sequences
        end
      end

      def dump(stream=$stdout)
        stream.puts "KyotoIndex dump: #{@path}"
        stream.puts
        db.cursor_process do |cur|
          stream.puts "== Metadata =="
          cur.jump('')
          while true
            k, v = cur.get(false)
            break if k[0] == "\xff"
            stream.puts "#{k}: #{v}"
            unless cur.step
              raise "could not advance cursor!"
            end
          end
          stream.puts "== Index records =="
          while pair = cur.get(true)
            _, chr, bin, s_start, s_end = pair[0].unpack(KEY_FMT)
            offset, len, text_size, n_seq, species_vec = pair[1].unpack(VAL_FMT)
            stream.puts "#{chr} [bin #{bin}] #{s_start}:#{s_end}"
            stream.puts "  offset #{offset}, length #{len}"
            stream.puts "  text size: #{text_size}"
            stream.puts "  sequences in block: #{n_seq}"
            stream.printf("  species vector: %016x\n", species_vec)
          end
        end
      end

      ## Retrieval:
      ##  1. merge the intervals of interest
      ##  2. for each interval, compute the bins with #bin_all
      ##  3. for each bin to search, make a list of intervals of
      ##     interest
      ##  4. compute the spanning interval for that bin
      ##  5. start at the beginning of the bin
      ##  6. if a record intersects the spanning interval: 
      ##    A. #find an interval it intersects
      ##    B. if found, add to the fetch list
      ##  7. if a record starts past the end of the spanning interval,
      ##     we are done scanning this bin.
      ##
      ## Optimizations:
      ##  * once we reach the start of the spanning interval,
      ##    all records start in it until we see a record starting
      ##    past it.
      ##  * as record starts pass the start of intervals of interest,
      ##    pull those intervals off the list

      # Build a fetch list of alignment blocks to read, given an array
      # of Bio::GenomicInterval objects
      def fetch_list(intervals, filter_spec={})
        to_fetch = []
        filter_spec ||= {}
        filters = Filters.build(filter_spec, self)
        chrom = intervals.first.chrom
        chrom_id = index_sequences[chrom]
        unless chrom_id
          raise "chromosome #{chrom} not indexed!"
        end
        if intervals.find { |i| i.chrom != chrom }
          raise "all intervals must be for the same chromosome!"
        end
        # for each bin, build a list of the intervals to look for there
        bin_intervals = Hash.new { |h, k| h[k] = [] }
        intervals.each do |i|
          i.bin_all.each { |bin| bin_intervals[bin] << i }
        end
        db.cursor_process do |cur|
          bin_intervals.each do |bin, bin_intervals_raw|
            bin_intervals = bin_intervals_raw.sort_by { |i| i.zero_start }
            # compute the start and end of all intervals of interest
            spanning_start = bin_intervals.first.zero_start
            spanning_end = bin_intervals.collect {|i| i.zero_end}.sort.last
            # scan from the start of the bin
            cur.jump(bin_start_prefix(chrom_id, bin))
            while pair = cur.get(true)
              @key.read(pair[0])
              @val.read(pair[1])
              if (@key.seq_id != chrom_id) \
                || (@key.bin != bin) \
                || @key.seq_start >= spanning_end
                # we've hit the next bin, or chromosome, or gone past
                # the spanning interval, so we're done with this bin
                break
              end
              if @key.seq_end >= spanning_start # possible overlap
                c_int = GenomicInterval.zero_based(chrom,
                                                   @key.seq_start,
                                                   @key.seq_end)
                if bin_intervals.find { |i| i.overlapped?(c_int) }
                  # TODO: simplify
                  pair_2 = [@key, @val]
                  if filters.match(pair_2)
                    to_fetch << [@val.b_offset, @val.b_length]
                  end
                end
              end
            end
          end # bin_intervals.each
        end # #cursor_process
        return to_fetch
      end # #fetch_list

     def build_default(parser)
        first_block = parser.parse_block
        ref_seq = first_block.sequences.first.source
        @index_sequences = { ref_seq => 0 }
        store_index_sequences!
        index_block(first_block)
        parser.each_block { |b| index_block(b) }
      end

      def load_index_sequences
        h = {}
        db.match_prefix("sequence:").each do |key|
          _, name = key.split(':', 2)
          id = db[key].to_i
          h[name] = id
        end
        @index_sequences = h
      end

      def store_index_sequences!
        index_sequences.each do |name, id|
          db.set("sequence:#{name}", id.to_s)
        end
      end

      def species_id_for_seq(seq)
        parts = seq.split('.')
        if parts.size == 2
          species_name = parts[0]
          if species.has_key? species_name
            return species[species_name]
          else
            species_id = @species_max_id + 1
            if species_id >= MAX_SPECIES
              raise "cannot index MAF file with more than #{MAX_SPECIES} species"
            end
            species[species_name] = species_id
            db["species:#{species_name}"] = species_id
            @species_max_id = species_id
            return species_id
          end
        else
          # not in species.sequence format, apparently
          return nil
        end
      end

      def index_block(block)
        entries_for(block).each do |k, v|
          db.set(k, v)
        end
      end

      def build_block_value(block)
        ## TODO: RSpec to ensure that all fields are reset
        @val.b_offset = block.offset
        @val.b_length = block.size
        @val.text_size = block.text_size
        @val.n_seq = block.sequences.size
        species_vec = BitString.new(0, 64)
        block.sequences.each do |seq|
          species_vec[species_id_for_seq(seq.source)] = 1
        end
        @val.species_vec = species_vec.to_i
      end

      def entries_for(block)
        ## TODO: RSpec to ensure that all fields are reset
        e = []
        build_block_value(block)
        block.sequences.each do |seq|
          seq_id = index_sequences[seq.source]
          next unless seq_id
          @key.seq_id = seq_id
          @key.seq_start = seq.start
          seq_end = seq.start + seq.size
          @key.seq_end = seq_end
          @key.bin = Bio::Ucsc::UcscBin.bin_from_range(seq.start, seq_end)
          e << [@key.to_binary_s, @val.to_binary_s]
        end
        return e
      end
    end # class KyotoIndex

    class Filter
      include KVHelpers

      def call(e)
        match(e)
      end
    end

    class AllSpeciesFilter < Filter
      attr_reader :bs
      def initialize(species, idx)
        bs = BitString.new(0, 64)
        species.each do |species_name|
          bs[idx.species.fetch(species_name)] = 1
        end
        @bs = bs
      end

      def match(entry)
        vec = extract_species_vec(entry)
        (@bs & vec) == @bs
      end
    end

    class AtLeastNSequencesFilter < Filter
      attr_reader :n
      def initialize(n, idx)
        @n = n
      end

      def match(entry)
        bs = BitString.new(extract_species_vec(entry))
        bs.population(1) >= n
      end
    end

    class MaxSizeFilter < Filter
      def initialize(n, idx)
        @n = n
      end
      def match(entry)
        extract_text_size(entry) <= @n
      end
    end

    class MinSizeFilter < Filter
      def initialize(n, idx)
        @n = n
      end
      def match(entry)
        extract_text_size(entry) >= @n
      end
    end

    class Filters
      include KVHelpers

      FILTER_CLASSES = {
        :with_all_species => MAF::AllSpeciesFilter,
        :at_least_n_sequences => MAF::AtLeastNSequencesFilter,
        :min_size => MAF::MinSizeFilter,
        :max_size => MAF::MaxSizeFilter
      }

      def self.build(spec, idx)
        l = spec.collect do |key, val|
          if FILTER_CLASSES.has_key? key
            FILTER_CLASSES[key].new(val, idx)
          else
            raise "Unsupported filter key #{key}!"
          end
        end
        return Filters.new(l)
      end

      def initialize(l)
        @l = l
      end

      def match(entry)
        return ! @l.find { |f| ! f.call(entry) }
      end
    end

  end # module MAF
  
end
