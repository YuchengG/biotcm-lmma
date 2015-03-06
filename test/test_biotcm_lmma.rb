require 'minitest/autorun'

$LOAD_PATH.unshift(File.expand_path('../../lib', __FILE__))
require 'biotcm-lmma'

describe BioTCM::Apps::LMMA do
  describe 'when used to construct LM-based network' do
    before do
      @lmma = BioTCM::Apps::LMMA.new
      @context = BioTCM::Databases::Medline.new('psoriasis') & 'breast cancer'
    end

    it 'could use co-occurrence information' do
      @lmma.lm(@context, method:'co-occurrence',
               node_savepath:'node.tab', edge_savepath:'edge.tab')
    end
  end

  describe 'when used to construct MA-based network' do
  end
end
