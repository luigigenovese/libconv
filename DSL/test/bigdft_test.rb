[ '../lib', 'lib' ].each { |d| $:.unshift(d) if File::directory?(d) }
require 'bigdft'
include BigDFT
require 'minitest/autorun'
require_relative 'bigdft/system'
require_relative 'bigdft/wave_function'
