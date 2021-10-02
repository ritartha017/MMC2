# frozen_string_literal: true

# Conjugate Gradient.
# Hestenes-Stiefel Algorithm
# Fletcher Reeves Algorithm
# Polak-Ribiere Algorithm

require 'rubocop'

@a = 4
@b = 3
@x = { x: @a, y: @b } # =>x0
@k = 0

def fetch_partial_derivatives
  [2 * @a * @x[:x] + 3 * @x[:y] - @a,
   3 * @x[:x] + 2 * @b * @x[:y] - @b]
end

def Q
  [2 * @a, 3,
   3, 2 * @b]
end

def fetch_gk
  fetch_partial_derivatives
end

def fetch_dk(g1, beta, d)
  [-g1[0] + beta * d[0],
   -g1[1] + beta * d[1]]
end

def fetch_next_x(alfa, d)
  @x[:x] = @x[:x] + alfa * d[0]
  @x[:y] = @x[:y] + alfa * d[1]
  @x
end

def gk_x_dk(g, d)
  g[0] * d[0] + g[1] * d[1]
end

def dk_Q_dk(_g, d)
  multiplix1 = [d[0] * Q()[0] + d[1] * Q()[1],
                d[0] * Q()[2] + d[1] * Q()[3]]
  multiplix1[0] * d[0] + multiplix1[1] * d[1]
end

def fletcher_reeves(g1, d)
  multiplix1 = [g1[0] * Q()[0] + g1[1] * Q()[1],
                g1[0] * Q()[2] + g1[1] * Q()[3]]
  numerator = multiplix1[0] * d[0] + multiplix1[1] * d[1]
  multiplix3 = [d[0] * Q()[0] + d[1] * Q()[1],
                d[0] * Q()[2] + d[1] * Q()[3]]
  denumerator = multiplix3[0] * d[0] + multiplix3[1] * d[1]
  numerator / denumerator.to_f
end

def hestenes_stiefel(g1, d)
  g1[0] * g1[0] + g1[1] * g1[1] / (d[0] * (-d[0]) + d[1] * (-d[1])).to_f
end

def polak_ribiere(g1, d)
  multiplix1 = [d[0] * Q()[0] + d[1] * Q()[1],
                d[0] * Q()[2] + d[1] * Q()[3]]
  denumerator = multiplix1[0] * d[0] + multiplix1[1] * d[1]
  g1[0] * g1[0] + g1[1] * g1[1] / denumerator.to_f
end

def minimize
  g = fetch_gk
  d = [-g[0], -g[1]]
  return if (g[0]).zero? && (g[1]).zero?

  until (g[0]).zero? && (g[1]).zero?
    @k += 1
    alfa = -(gk_x_dk(g, d) / dk_Q_dk(g, d).to_f)
    @x = fetch_next_x(alfa, d)
    g = fetch_gk
    beta = polak_ribiere(g, d)
    d = fetch_dk(g, beta, d)
  end
  @x
end

puts "\nx* = #{minimize}"
puts "N iterations = #{@k}"
