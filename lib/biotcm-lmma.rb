require 'biotcm'
require 'optparse'

# A combined literature mining and microarray analysis approach to construct
# gene networks of a specific biological system.
class BioTCM::Apps::LMMA < BioTCM::Apps::App
  include BioTCM::Modules::WorkingDir

  # Current version
  VERSION = "0.1.0"

  # Terminal entry
  def run
    OptionParser.new do |opts|
      opts.banner = 'Usage: biotcm lmma [OPTIONS]'
    end.parse!

    fail NotImplementedError
  end
  # Build LM-based network
  # @param context [BioTCM::Databases::Medline]
  #   a medline object representing the context for LM-based network construction
  # @param opts [Hash]
  # @option opts [String] :method only 'co-occurrence' by now
  # @option opts [String] :node_savepath
  # @option opts [String] :edge_savepath
  # @return [Hash]
  def lm(context, opts = {})
    raise ArgumentError unless context.is_a?(BioTCM::Databases::Medline)
    lm_cooccurrence(context, opts)
    return self
  end

  private

  # Use gene co-occurrences to build LM-based network
  def lm_cooccurrence(context, opts = {})
    path_to_medlines =
      BioTCM.path_to("tmp/lmma_medlines_#{BioTCM.get :stamp}.txt")

    # Download abstracts
    context.download_abstracts(path_to_medlines)

    @lmnet_node = {}
    @lmnet_edge = {}
    counter  = 0
    pmid     = nil
    abstract = ''
    detector = BioTCM::Apps::GeneDetector.new

    f_abstracts = File.open(path_to_medlines)
    f_abstracts.each do |line|
      if line =~ /^PMID- +(\d+)/
        pmid = $1
        counter += 1
        BioTCM.logger.info('LMMA') { "Analyzing article \##{counter}/#{context.count}..." }
      elsif line =~ /^AB  -/
        abstract = line.gsub(/^AB  -\s*/, '').chomp
        abstract += line.gsub(/^\s+/, ' ') while (line = f_abstracts.gets.chomp) =~ /^\s/

        # Split into sentences
        sentences = abstract.split(/[.?!][\s]?/)

        # Identify genes
        gene_sets = sentences.collect { |sentence| detector.detect(sentence) }

        # Update lmnet nodes
        gene_sets.flatten.uniq.each do |gene|
          if @lmnet_node[gene]
            @lmnet_node[gene] << pmid
          else
            @lmnet_node[gene] = [pmid]
          end
        end

        # Update lmnet edges
        gene_sets.collect do |gene_set|
          gene_set.product(gene_set)
        end.inject([]) {|a, b| a + b}.collect do |gene_pair|
          gene_pair.sort
        end.uniq.each do |gene_pair|
          next if gene_pair[0] == gene_pair[1]
          key = gene_pair.join("\t")
          if @lmnet_edge[key]
            @lmnet_edge[key] << pmid
          else
            @lmnet_edge[key] = [pmid]
          end
        end
      end
    end

    # Export
    if(opts[:node_savepath])
      tab = Table.new(primary_key: 'Symbol',
                      col_keys:   ['Occurrence count', 'PMID(s)'])
      @lmnet_node.keys.sort_by { |k| - @lmnet_node[k].size }.each do |k|
        tab.row(k, [@lmnet_node[k].size, @lmnet_node[k].join(', ')])
      end
      tab.save(opts[:node_savepath])
    end
    if(opts[:edge_savepath])
      tab = Table.new(primary_key: "Source\tTarget",
                      col_keys:   ['Co-occurrence count', 'PMID(s)'])
      @lmnet_edge.keys.sort_by { |k| - @lmnet_edge[k].size }.each do |k|
        tab.row(k, [@lmnet_edge[k].size, @lmnet_edge[k].join(', ')])
      end
      tab.save(opts[:edge_savepath])
    end
  end
end
