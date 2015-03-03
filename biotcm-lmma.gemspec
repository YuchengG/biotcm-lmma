$:.unshift File.expand_path("../lib", __FILE__)
require 'rake'
require 'biotcm-lmma'

Gem::Specification.new do |s|
  s.name        = "biotcm-lmma"
  s.platform    = Gem::Platform::RUBY
  s.summary     = "A combined literature mining and microarray analysis approach."
  s.description = "A combined literature mining and microarray analysis approach to construct gene networks of a specific biological system."

  s.version     = BioTCM::Apps::LMMA::VERSION
  s.license     = 'MIT'

  s.authors     = ["Aidi Stan"]
  s.email       = ["aidistan@live.cn"]
  s.homepage    = 'https://github.com/biotcm/biotcm-lmma'

  s.files       = FileList['lib/**/*'].to_a

  s.required_ruby_version = '>= 2.0.0'
  s.add_runtime_dependency "biotcm", "~> 0.2.0"
end
