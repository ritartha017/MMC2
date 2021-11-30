# frozen_string_literal: true

require 'rulp'

Rulp.log_level = Logger::DEBUG

variables = [Oatmeal_i, Chicken_i, Eggs_i, Whole_milk_i, Cherry_pie_i, Pork_and_beans_i]

given = Array[variables.map { |x| x >= 0 }]

constraints = [
  110 * Oatmeal_i + 205 * Chicken_i + 160 * Eggs_i + 160 *
                                                     Whole_milk_i + 420 * Cherry_pie_i + 260 * Pork_and_beans_i >= 2000,
  4 * Oatmeal_i + 32 * Chicken_i  + 13 * Eggs_i + 8 *
                                                  Whole_milk_i + 4 * Cherry_pie_i + 14 * Pork_and_beans_i >= 55,
  2 * Oatmeal_i + 12 * Chicken_i  + 54 * Eggs_i + 285 *
                                                  Whole_milk_i + 22 * Cherry_pie_i + 80 * Pork_and_beans_i >= 800
]

objective = 0.3 * Oatmeal_i + 2.40 * Chicken_i + 1.30 * Eggs_i + 0.90 * Whole_milk_i + 2.0 * Cherry_pie_i + 1.9 * Pork_and_beans_i

Rulp::Min(objective)[constraints].glpk

result = objective.evaluate

print "\nMin cost = $#{result}"
