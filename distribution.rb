# encoding: utf-8
require 'mathn.rb'

module DistributionTask

  def distribution_matrix x, f
    b = [f.first.zip(x)] 
    (1...f.length).each do |i|
      b << []
      x.each do |y|
	it = x.select{|item| item <= y}.map{|z|
	  [b[i-1][x.index(y-z)].first + f[i][x.index(z)], z] }.group_by{|item| item.first}.max.last
	b[i] << ([it.first.first, it.map{|item| item[1]}])
      end
    end
    b
  end
  
  def print_matrix matrix = []
    puts matrix.map{|item| item.map{|i| "#{i.first} / #{i[1]}"}.join(" | ")}.join("\n")
  end
  
end

include DistributionTask

#x = [0,10,20,30,40,50,60]

#f = [
#  [0,15,20,25,30,35,40],
#  [0,10,15,19,25,28,33],
#  [0,10,18,20,25,29,35],
#  [0,15,20,25,30,35,40]
#  ]

x = [0,5,10,15,20,25, 30]

f = [
  [0,10,20,29,37,45,50],
  [0,15,24,30,39,45,49],
  [0,14,23,32,41,50,53]
  ]

#f = [
#  [0, 18, 36, 54, 72, 90],
#  [0, 20, 40, 60, 80, 89],
#  [0, 22, 45, 66, 80, 92],
#  [0, 15, 30, 45, 60, 75]
#  ]
print_matrix distribution_matrix x, f