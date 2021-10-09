# frozen_string_literal: true

# Conjugate Gradient.
# Hestenes-Stiefel Algorithm
# Fletcher Reeves Algorithm
# Polak-Ribiere Algorithm

# require 'rubocop'

@a = 4
@b = 3

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

  res = gk_dk / (dk_q[0] * dk[0] + dk_q[1] * dk[1]).to_f
  res.nan? || res.infinite? ? 0 : res
end

def fetch_x(xk, alphak, dk)
  [xk[0] + alphak * dk[0],
   xk[1] + alphak * dk[1]]
end

def fetch_beta(gk1, dk)
  gk1_q = [gk1[0] * Q[0] + gk1[1] * Q[1],
           gk1[0] * Q[2] + gk1[1] * Q[3]]
  gk1_q_dk = gk1_q[0] * dk[0] + gk1_q[1] * dk[1]
  dk_q = [dk[0] * Q[0] + dk[1] * Q[1], dk[0] * Q[2] + dk[1] * Q[3]]
  dk_q_dk = dk[0] * dk_q[0] + dk[1] * dk_q[1]

  res = gk1_q_dk / dk_q_dk.to_f
  res.nan? || res.infinite? ? 0 : res
end

def fetch_dk1(gk1, beta, dk)
  [-gk1[0] + beta * dk[0],
   -gk1[1] + beta * dk[1]]
end

def minimize
  alphak = xk1 = beta = k = 0
  xk = [@a, @b]
  dk1 = []
  gk1 = fetch_partial_derivatives(xk[0], xk[1])
  dk1[0] = -gk1[0]
  dk1[1] = -gk1[1]
  p "G= #{gk1}"
  loop do
    k += 1
    alphak = -fetch_alpha(gk1, dk1)
    xk1 = fetch_x(xk, alphak, dk1)
    gk1 = fetch_partial_derivatives(xk1[0], xk1[1])
    beta = fetch_beta(gk1, dk1)
    puts    "k= #{k} \n"    + "alpha= #{alphak}\n" \
            "x= #{xk1} \n"  + "G= #{gk1} \n"       \
            "B= #{beta}"    + "\n\n"
    dk1 = fetch_dk1(gk1, beta, dk1)
    xk = xk1
    return xk if (gk1[0]).zero? && (gk1[1]).zero?
  end
end

puts "\nx* = #{minimize}"
