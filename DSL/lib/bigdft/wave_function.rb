module BOAST

  class WaveFunction
    attr_accessor :data_space
    attr_reader   :system
    attr_reader   :spaces
    attr_reader   :leading_dimensions

    def initialize(system, **options)
      @system = system
      @space = options.fecth(:spaces, [@system.reference_space]*@system.dimension)
    end

  end

end
