require 'biotcm'
require 'optparse'

# LMMA
class BioTCM::Apps::LMMA < BioTCM::Apps::App
  include BioTCM::Modules::WorkingDir

  # Version
  VERSION = "0.0.1"

  # Terminal entry
  def run
    OptionParser.new do |opts|
      opts.banner = "Usage: biotcm lmma [OPTIONS]"
    end.parse!

    fail NotImplementedError
  end
end
