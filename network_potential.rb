# encoding: utf-8
require 'matrix.rb'
require 'mathn.rb'

module PotentialMethodHelper

  def find_path trade, network, cost, x_init = nil, basis_init = nil
    x_init, basis_init  = *find_init_basis(trade, network) if x_init.nil? || basis_init.nil?
    puts "Second stage: "
    print "Cost: ", cost
    x = deep_copy x_init
    print "X: ", x
    basis = deep_copy basis_init
    print "Basis: ", basis, false
    optimal = false
    counter = 1
    until optimal
      puts "Iteration ##{counter}"
      counter += 1
      x, basis, optimal = *iterate(basis, cost, network, x)
    end
    print "Optimal x: ", x
    print "Optimal basis: ", basis, false
    [x, basis]
  end
  
  #private
 
  def iterate basis, cost, network, x
    potentials = count_potentials basis, cost
    puts potentials.each_with_index.map{|u, i| "u#{i+1} = #{u}"}
    estimations = count_estimations basis, network, potentials, cost
    puts "Check optimality criterion..."
    return [x, basis, true] if optimal? estimations
    puts "Optimality criterion is not satisfied"
    index0 = choose_not_optimal_indexes basis, estimations
    puts "(i0, j0) = (#{index0[0]+1}, #{index0[1]+1})"
    step0, indexz = *make_step0(index0, basis, x)
    raise "No solution" if indexz.nil?
    puts "step0 = #{step0}, (i*, j*) = (#{indexz[0]+1}, #{indexz[1]+1})"
    x_new = make_x_new x, step0, index0, basis
    print "New x: ", x_new
    basis_new = make_basis_new basis, index0, indexz
    print "New basis: ", basis_new, false
    [x_new, basis_new, false]
  end
  
  def make_basis_new basis, index0, indexz
    basis_new = deep_copy basis
    basis_new[index0[0]][index0[1]] = 1
    basis_new[indexz[0]][indexz[1]] = 0
    basis_new
  end
  
  def print phrase, array, print_value = true
    puts phrase
    array.each_with_index{|line, i| line.each_with_index{|item, j| puts "#{i+1} -> #{j+1} #{": #{item}" if print_value}" unless item==0}}
  end
  
  def make_x_new x, step0, index0, basis
    x_new = deep_copy x
    parents = find_cycle basis, index0
    i = index0[0]
    until i==index0[1]
      #puts "#{i+1} #{parents[i]+1}"
      if basis[i][parents[i]]==1
	x_new[i][parents[i]] -= step0
      else
	x_new[parents[i]][i] += step0
      end
      i = parents[i]
    end
    x_new[index0[0]][index0[1]] += step0
    x_new
  end
  
  def deep_copy obj
    Marshal.load(Marshal.dump obj)
  end
  
  def make_step0 index0, basis, x
    puts "Count steps for basis edges and choose step0: "
    step0 = [Float::INFINITY]
    parents = find_cycle basis, index0
    #p parents
    i = index0[0]
    until i==index0[1]
      step0 = [x[i][parents[i]],[i,parents[i]]] if step0[0] > x[i][parents[i]] && basis[i][parents[i]]==1
      puts "step#{i+1}#{parents[i]+1} = #{basis[i][parents[i]]==1 ? x[i][parents[i]] : Float::INFINITY}"
      i = parents[i]
    end
    step0
  end
  
  def find_cycle basis, index0
    basis_ext = deep_copy basis
    basis_ext = make_not_oriented basis_ext
    basis_ext.size.times{|i| basis_ext[index0[0]][i] = 0}
    basis_ext[index0[0]][index0[1]] = 1
    visited = [false]*basis_ext.size
    parents = [-1]*basis_ext.size
    find_deep index0[0], visited, basis_ext, parents
    parents
  end
  
  def make_not_oriented matrix
    m = deep_copy matrix
    m.each_with_index{|line, i| line.each_with_index{|item, j| m[j][i] = 1 if item==1}}
    m
  end
  
  def find_deep point, visited, basis, parents
    visited[point] = true
    basis[point].each_with_index do |item, i|
      to = basis[point][i]
      #puts "#{point+1} -> #{i+1}, #{visited[i]}, #{to}"
      if !visited[i] && to==1
	parents[i] = point
	return true if find_deep i, visited, basis, parents
      elsif parents[point]!=i && to==1
	parents[i] = point
	return true
      end
    end
    false
  end
  
  def optimal? estimations
    estimations.each_with_index.none? do 
      |line, i| line.each_with_index.any? do |estimate, j|
	puts "Not optimal delta#{i+1}#{j+1} = #{estimate}" if !estimate.nil? && estimate < 0
	estimate.nil? ? false : estimate < 0
      end
    end
  end
  
  def choose_not_optimal_indexes basis, estimations
    puts "Choose i0, j0 with the lowest delta:"
    min = []
    estimations.each_with_index do |line, i| 
      line.each_with_index do |edge, j| 
	min = [edge,[i,j]]if basis[i][j]==0 && !edge.nil? && (min.empty? || edge < min[0] ) && edge < 0
      end
    end
    min[1]
  end
  
  def count_potentials basis, cost
    puts "Count potentials for basis edges: "
    matrix = []
    b = []
    basis.each_with_index do |line, i|
      line.each_with_index do |c, j| 
	if c == 1
	  puts "u#{i+1} - u#{j+1} = c#{i+1}#{j+1} = #{cost[i][j]}"
	  line = [0]*cost.size
	  line[i], line[j] = 1, -1
	  matrix <<  line
	  b << cost[i][j]
	end
      end
    end
    potentials = [0]+(Matrix.rows(make_first_null(matrix)).inverse*Vector.elements(b)).to_a 
  end
  
  def count_estimations basis, network, potentials, cost
    puts "Count estimations for not basis edges: "
    puts "delta_ij = c_ij - (u_i -u_j), (i,j) not basis"
    estimations = Array.new(basis.size){[nil]*basis.size}
    network.each_with_index do |line, i|
      line.each_with_index do |edge, j|
	estimations[i][j] = cost[i][j] - (potentials[i] - potentials[j])  if edge!=0 && basis[i][j]==0 
	puts "delta#{i+1}#{j+1} = #{cost[i][j]} - (#{potentials[i]} - #{potentials[j]}) = #{estimations[i][j]}" if edge!=0 && basis[i][j]==0
      end
    end
    estimations
  end
  
  def make_first_null matrix
    matrix.map{|line| line[1,line.size]}
  end
  
  def find_init_basis trade, network
    puts "First stage: "
    network_init = make_network_init trade, network
    print "Initial network: ", network_init, false
    cost_init = make_cost_init trade
    print "Initial cost: ", cost_init
    basis_init = deep_copy cost_init
    print "Initial basis: ", basis_init, false
    x_init = make_x_init trade
    print "Initial x: ", x_init
    x = deep_copy x_init
    basis = deep_copy basis_init
    cost = deep_copy cost_init
    optimal = false
    counter = 1
    until optimal
      puts "Iteration ##{counter}"
      counter += 1
      x, basis, optimal = *iterate(basis, cost, network_init, x)
    end
    print "Optimal x: ", x
    print "Optimal basis: ", basis, false
    x, basis = *delete_artificials(x, basis)
    [x, basis]
  end
  
  def delete_artificials x, basis
    x = deep_copy(x)
    x.pop
    basis = deep_copy(basis)
    basis.pop
    [x.map{|line| line[0, line.size-1]},basis.map{|line| line[0, line.size-1]} ]
  end
  
  def make_x_init trade
    x_init =  Array.new(trade.size+1){[0]*(trade.size+1)}
    trade.each_with_index do |tr, i|
      if tr > 0
	x_init[i][trade.size] = tr
      elsif
	x_init[trade.size][i] = tr.abs
      end
    end
    x_init
  end
  
  def make_cost_init trade
    trade.map do |cost|
      [0]*trade.size << (cost>= 0 ? 1 : 0)
    end << (trade.map{|cost| cost < 0 ? 1 : 0}+[0])
  end
  
  def make_network_init trade, network
    network.each_with_index.map do |point, i|
      point << (trade[i] >= 0 ? 1 : 0)
    end << (trade.map{|cost| cost < 0 ? 1 : 0} << 0)
  end
  
end

include PotentialMethodHelper

network = [
  [0,0,0,0,1,1,1,1,1],
  [0,0,0,0,1,1,1,1,1],
  [0,0,0,0,1,1,1,1,1],
  [0,0,0,0,1,1,1,1,1],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0]
  ]
cost = [
  [0,0,0,0,8,10,8,10,4],
  [0,0,0,0,3,10,1,6,1],
  [0,0,0,0,10,9,10,8,3],
  [0,0,0,0,1,12,10,1,10],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0]
  ]
basis_init = [
  [0,0,0,0,0,1,0,1,1],
  [0,0,0,0,0,0,1,0,1],
  [0,0,0,0,0,0,0,0,1],
  [0,0,0,0,1,0,0,1,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0]
  ]
x_init = [
  [0,0,0,0,0,6,0,4,0],
  [0,0,0,0,0,0,8,0,3],
  [0,0,0,0,0,0,0,0,8],
  [0,0,0,0,6,0,0,3,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0,0]
  ]
trade = [10,11,8,9,-6,-6,-8,-7,-11]

find_path trade, network, cost, x_init, basis_init
