require 'minitest/autorun'
require 'tempfile'

$LOAD_PATH.unshift(File.expand_path('../../lib', __FILE__))
require 'biotcm-lmma'

# Prepare data
lm_data = Tempfile.new('medline')
medline = BioTCM::Databases::Medline.new('psoriasis') & 'breast cancer'
medline.download_abstracts(lm_data.path)
ma_datapath = File.expand_path('../test_biotcm_lmma/fake.tab', __FILE__)

describe BioTCM::Apps::LMMA do
  before do
    @lmma = BioTCM::Apps::LMMA.new
  end

  it 'could build LM-based network' do
    assert(@lmma.lm(lm_data.path, total_number: medline.count))
    assert_instance_of(Table, @lmma.lm_node)
    assert_instance_of(Table, @lmma.lm_edge)
  end

  it 'could build MA-based network' do
    assert(@lmma.ma(ma_datapath))
    assert_instance_of(Table, @lmma.ma_edge)
  end

  it 'could build LMMA-based network' do
    assert(@lmma.lmma(lm_data.path, ma_datapath, total_number: medline.count))
    assert_instance_of(Table, @lmma.lmma_edge)
  end
end
