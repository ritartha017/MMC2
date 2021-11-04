# frozen_string_literal: true

# Gradient Descent Algorithm applied in Ruby.

require 'rubocop'
require 'yaml'

EPSELON = 10**-6
SIGMA = 0.06
GAMMA = 0.5

@a = 4
@b = 3
@x = { x: @a, y: @b } # =>x0
@lamda = 1
@k = 0

def fetch_partial_derivatives(x, y)
  { x: 2 * @a * x + 3 * y - @a,
    y: 3 * x + 2 * @b * y - @b }
end

def substitute_into_fxy(x, y)
  @a * x**2 + 3 * x * y + @b * (y**2) - @a * x - @b * y
end

def next_point(xk, delta_fx)
  { x: xk[:x] - @lamda * (delta_fx[:x]),
    y: xk[:y] - @lamda * (delta_fx[:y]) }
end

def get_magnitude(x, y)
  Math.sqrt(x**2 + y**2)
end

def pp_print
  out = { "x#{@k}" => [{ xcoordinate: @x, fx: @fx, Δx: @delta_fx, magnitude: @magnitude_fx },
                       { zcoordinate: @z, fz: @fz, Δy: @delta_fz, magnitude: @magnitude_fz },
                       lamda: @lamda, condition: @condition] }
  puts out.to_yaml
end

def fetch_local_minimum
  loop do
    @fx = substitute_into_fxy(@x[:x], @x[:y])
    @delta_fx = fetch_partial_derivatives(@x[:x], @x[:y])
    @magnitude_fx = get_magnitude(@delta_fx[:x], @delta_fx[:y])

    loop do
      @z = next_point(@x, @delta_fx)
      @fz = substitute_into_fxy(@z[:x], @z[:y])
      @delta_fz = fetch_partial_derivatives(@z[:x], @z[:y])
      @magnitude_fz = get_magnitude(@delta_fz[:x], @delta_fz[:y])

      @condition = @fz - @fx <= -SIGMA * @lamda * @magnitude_fx**2
      @lamda *= GAMMA unless @condition
      pp_print
      break unless @condition != true
    end

    @x = @z
    @k += 1
    break if @magnitude_fz < EPSELON
  end
  @fz
end

stationary_point = fetch_local_minimum
