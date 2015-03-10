require 'minitest/autorun'

$LOAD_PATH.unshift(File.expand_path('../../lib', __FILE__))
require 'biotcm-lmma'

describe BioTCM::Apps::LMMA do
  before do
    @lmma     = BioTCM::Apps::LMMA.new
    @context  = BioTCM::Databases::Medline.new('psoriasis') & 'breast cancer'
    @datapath = File.expand_path('../test_biotcm_lmma/fake.tab', __FILE__)
  end

  it 'could build LM-based network' do
    assert(@lmma.lm(@context))
    assert_instance_of(Table, @lmma.lm_node)
    assert_instance_of(Table, @lmma.lm_edge)
  end

  it 'could build MA-based network' do
    assert(@lmma.ma(@datapath))
    assert_instance_of(Table, @lmma.ma_edge)
  end

  it 'could build LMMA-based network' do
    assert(@lmma.lmma(@context, @datapath))
    assert_instance_of(Table, @lmma.lmma_edge)
  end
end
