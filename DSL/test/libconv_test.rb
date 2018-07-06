[ '../lib', 'lib' ].each { |d| $:.unshift(d) if File::directory?(d) }
require 'libconv'

require 'libconv/space'
require 'libconv/system'
require 'libconv/wave_function'
include LibConv

require 'minitest/autorun'
require_relative 'libconv/system'
require_relative 'libconv/wave_function'

