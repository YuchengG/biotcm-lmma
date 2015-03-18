require 'erb'
require 'biotcm'
require 'optparse'
require 'tempfile'

# A combined literature mining and microarray analysis approach to construct
# gene networks of a specific biological system.
class BioTCM::Apps::LMMA < BioTCM::Apps::App
  include BioTCM::Modules::WorkingDir

  # Current version
  VERSION = "0.2.0"

  attr_reader :lm_node, :lm_edge
  attr_reader :ma_edge
  attr_reader :lmma_edge

  # Terminal entry
  def run
    OptionParser.new do |opts|
      opts.banner = 'Usage: biotcm lmma [OPTIONS]'
    end.parse!

    fail NotImplementedError
  end
  # Build LM-based network
  # @param datapath [String] a path to a downloaded medline file
  # @param opts [Hash]
  # @option opts [String] :node_savepath a path to save node table
  # @option opts [String] :edge_savepath a path to save edge table
  # @option opts [Fixnum] :total_number the total number of the articles to analyze
  # @return [self]
  def lm(datapath, opts = {})
    raise 'File not exists' unless File.exist?(datapath)

    lm_cooccurrence(datapath, opts)

    return self
  end
  # Build MA-based network
  # @param datapath [String] a path to a gene expression file
  # @param opts [Hash]
  # @return [self]
  def ma(datapath, opts = {})
    raise 'File not exists' unless File.exist?(datapath)

    ma_pairwise(datapath, opts)

    return self
  end
  # Build LMMA-based network
  # @see #lm
  # @see #ma
  def lmma(lm_datapath, ma_datapath, lm_opts = {}, ma_opts = {})
    lm(lm_datapath, lm_opts)

    # Save LM-based network as the background network
    f_bgnet = Tempfile.new('bgnet')
    f_bgnet.write(@lm_edge.to_s)
    f_bgnet.rewind
    ma_opts[:bgnet_path] = f_bgnet.path

    ma(ma_datapath, ma_opts)

    # Build LMMA-based network
    @lmma_edge = Table.new(
      primary_key: "Source\tTarget",
      col_keys: @ma_edge.col_keys + @lm_edge.col_keys
    )
    @ma_edge.row_keys.each do |rkey|
      @lmma_edge.row(rkey, @ma_edge.row(rkey).merge(@lm_edge.row(rkey)))
    end

    return self
  end

  private

  # Use gene co-occurrences to build LM-based network
  def lm_cooccurrence(datapath, opts = {})
    node = {}
    edge = {}
    counter  = 0
    pmid     = nil
    abstract = ''
    detector = BioTCM::Apps::GeneDetector.new

    f_abstracts = File.open(datapath, 'r:utf-8')
    f_abstracts.each do |line|
      if line =~ /^PMID- +(\d+)/
        pmid = $1
        counter += 1
        BioTCM.logger.info('LMMA') { "Analyzing article \##{counter}#{
          opts[:total_number] ? "/#{opts[:total_number]}" : ''
        }..." }
      elsif line =~ /^AB  -/
        abstract = line.gsub(/^AB  -\s*/, '').chomp
        abstract += line.gsub(/^\s+/, ' ') while (line = f_abstracts.gets.chomp) =~ /^\s/

        # Split into sentences
        sentences = abstract.split(/[.?!][\s]?/)

        # Identify genes
        gene_sets = sentences.collect { |sentence| detector.detect(sentence) }

        # Update nodes
        gene_sets.flatten.uniq.each do |gene|
          if node[gene]
            node[gene] << pmid
          else
            node[gene] = [pmid]
          end
        end

        # Update edges
        gene_sets.collect do |gene_set|
          gene_set.product(gene_set)
        end.inject([]) {|a, b| a + b}.collect do |gene_pair|
          gene_pair.sort
        end.uniq.each do |gene_pair|
          next if gene_pair[0] == gene_pair[1]
          key = gene_pair.join("\t")
          if edge[key]
            edge[key] << pmid
          else
            edge[key] = [pmid]
          end
        end
      end
    end

    # Build tables
    @lm_node = Table.new(
      primary_key: 'Symbol',
      col_keys:   ['Occurrence count', 'PMID(s)']
    )
    node.keys.sort_by { |k| - node[k].size }.each do |k|
      @lm_node.row(k, [node[k].size, node[k].join(', ')])
    end
    @lm_edge = Table.new(
      primary_key: "Source\tTarget",
      col_keys:   ['Co-occurrence count', 'PMID(s)']
    )
    edge.keys.sort_by { |k| - edge[k].size }.each do |k|
      @lm_edge.row(k, [edge[k].size, edge[k].join(', ')])
    end

    # Export
    @lm_node.save(opts[:node_savepath]) if(opts[:node_savepath])
    @lm_edge.save(opts[:edge_savepath]) if(opts[:edge_savepath])
  end
  # Use pairwise method to build MA-based network
  def ma_pairwise(datapath, opts = {})
    f_script = Tempfile.new('script')
    f_result = Tempfile.new('result')

    tp = File.read(File.expand_path('../biotcm-lmma/pairwise.R.erb', __FILE__))
    f_script.write(ERB.new(tp).result(binding))
    f_script.rewind
    `Rscript #{f_script.path}`

    @ma_edge = Table.new(primary_key: "Source\tTarget")
    colname = f_result.gets.chomp.split("\t")[2]
    f_result.each do |line|
      col = line.chomp.split("\t")
      next if @ma_edge.row("#{col[0]}\t#{col[1]}")
      next if @ma_edge.row("#{col[1]}\t#{col[0]}")
      @ma_edge.ele(col[0..1].sort.join("\t"), colname, col[2])
    end
  end
end
