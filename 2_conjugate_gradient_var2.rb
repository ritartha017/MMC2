# frozen_string_literal: true

@a = 4
@b = 3

def substitute_into_fxy(x, y)
  @a * x**2 + 3 * x * y + @b * (y**2) - @a * x - @b * y
end

def fetch_partial_derivatives(x, y)
  [2 * @a * x + 3 * y - @a,
   3 * x + 2 * @b * y - @b]
end

Q = [2 * @a, 3,
     3, 2 * @b].freeze

def fetch_alpha(gk, dk)
  gk_dk = gk[0] * dk[0] + gk[1] * dk[1]
  dk_q = [
    dk[0] * Q[0] + dk[1] * Q[1],
    dk[0] * Q[2] + dk[1] * Q[3]
  ]
  res = gk_dk / (dk_q[0] * dk[0] + dk_q[1] * dk[1]).to_f.send(:*, -1)
end

def fetch_x(xk, alphak, dk)
  [xk[0] + alphak * dk[0],
   xk[1] + alphak * dk[1]]
end

def get_hestenes_stiefel_beta(magnitude_k, magnitude_k1)
  res = (magnitude_k1**2).abs / (magnitude_k**2).abs.to_f
end

def get_fletcher_reeves_beta(xk, gk1)
  gk = fetch_partial_derivatives(xk[0], xk[1])
  (gk1[0] * gk1[0] + gk1[1] * gk1[1]) / (gk[0] * gk[0] + gk[1] * gk[1]).to_f
end

def fetch_dk1(gk1, beta, dk)
  [-gk1[0] + beta * dk[0],
   -gk1[1] + beta * dk[1]]
end

def get_magnitude(x, y)
  Math.sqrt(x**2 + y**2)
end

def minimize(a, b)
  alphak = xk1 = beta = k = 0
  xk = [a, b]
  gk1 = fetch_partial_derivatives(xk[0], xk[1])
  dk1 = [-gk1[0], -gk1[1]]
  magnitude_k = get_magnitude(gk1[0], gk1[1])
  puts "g0= #{gk1}\n" + "mk0= #{magnitude_k}"
  loop do
    k += 1
    alphak = fetch_alpha(gk1, dk1)
    xk1 = fetch_x(xk, alphak, dk1)
    fxk = substitute_into_fxy(xk1[0], xk1[1])
    gk1 = fetch_partial_derivatives(xk1[0], xk1[1])
    magnitude_k1 = get_magnitude(gk1[0], gk1[1])
    # beta = get_hestenes_stiefel_beta(magnitude_k, magnitude_k1)
    beta = get_fletcher_reeves_beta(xk, gk1)
    dk1 = fetch_dk1(gk1, beta, dk1)
    puts "k= #{k} \n" + "alpha= #{alphak}\n" + "x = #{xk1} \n" + "g = #{gk1} \n" \
         "B = #{beta}\n" + "Fx= #{fxk}" + "\nM = #{magnitude_k1}" + "\nd = #{dk1}\n\n"
    xk = xk1
    magnitude_k = magnitude_k1
    return xk if (gk1[0]).zero? && (gk1[1]).zero?
  end
end

puts "\nx* = #{minimize(@a, @b)}"
