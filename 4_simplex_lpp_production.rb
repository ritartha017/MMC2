# frozen_string_literal: true

require 'rulp'

Rulp.log_level = Logger::DEBUG

variables = [TN_i, VA_i, IPS_i, OLED_i]

given = Array[variables.map { |x| x >= 0 }]

constraints = [
  TN_i + VA_i + IPS_i + OLED_i <= 16,
  6 * TN_i + 5 * VA_i + 4 * IPS_i + 3 * OLED_i <= 110,
  4 * TN_i + 6 * VA_i + 10 * IPS_i + 13 * OLED_i <= 100
]

objective = 60 * TN_i + 70 * VA_i + 120 * IPS_i + 130 * OLED_i

Rulp::Max(objective)[constraints].glpk

result = objective.evaluate

print "\nMax profit = $#{result}"
