require 'pp'
require_relative 'BigDFT'
include BigDFT

psi = WaveFunction::new([32, 42, 28], :s1)


pp psi.to_real.to_wavelet
