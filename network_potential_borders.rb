# encoding: utf-8
require 'matrix.rb'
require 'mathn.rb'

module PotentialMethodHelper

  def find_path trade, network, cost, borders, x_init = nil, basis_init = nil
    x_init, basis_init  = *find_init_basis(trade, network, borders) if x_init.nil? || basis_init.nil?
    puts "Second stage: "
    x = deep_copy x_init
    basis = deep_copy basis_init
    optimal = false
    counter = 1
    until optimal
      puts "Iteration ##{counter}"
      counter += 1
      x, basis, optimal = *iterate(basis, cost, network, x, borders)
    end
    print "Optimal x: ", x
    print "Optimal basis: ", basis, false
    [x, basis]
  end
  
  #private
 
  def iterate basis, cost, network, x, borders
    potentials = count_potentials basis, cost
    puts potentials.each_with_index.map{|u, i| "u#{i+1} = #{u}"}
    estimations = count_estimations basis, network, potentials, cost
    puts "Check optimality criterion..."
    return [x, basis, true] if optimal? estimations, x, borders
    puts "Optimality criterion is not satisfied"
    index0 = choose_not_optimal_indexes basis, estimations, x, borders
    puts "(i0, j0) = (#{index0[0]+1}, #{index0[1]+1})"
    step0, indexz = *make_step0(index0, basis, x, borders)
    raise "No solution" if indexz.nil?
    puts "step0 = #{step0}, (i*, j*) = (#{indexz[0]+1}, #{indexz[1]+1})"
    x_new = make_x_new x, step0, index0, basis, borders
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
  
  def deep_copy obj
    Marshal.load(Marshal.dump obj)
  end
  
  def make_x_new x, step0, index0, basis, borders
    x_new = deep_copy x
    parents = find_cycle basis, index0
    direction_straight = x[index0[0]][index0[1]] == borders[index0[0]][index0[1]][0] ? true : false
    i = index0[0]
    until i==index0[1]
      if direction_straight
	if basis[i][parents[i]]==1 #обратная дуга
	  x_new[i][parents[i]] -= step0
	elsif basis[parents[i]][i]==1
	  x_new[parents[i]][i] += step0
	end
      else
	if basis[i][parents[i]]==1
	  x_new[i][parents[i]] += step0
	elsif basis[parents[i]][i]==1
	  x_new[parents[i]][i] -= step0
	end
      end                
      i = parents[i]
    end
    step = nil
    if direction_straight
      x_new[index0[0]][index0[1]] += step0
    else
      x_new[index0[0]][index0[1]] -= step0
    end
    x_new
  end

  
  def make_step0 index0, basis, x, borders
    puts "Count steps for basis edges and choose step0: "
    step0 = [Float::INFINITY]
    parents = find_cycle basis, index0
    direction_straight = x[index0[0]][index0[1]] == borders[index0[0]][index0[1]][0] ? true : false
    puts "Direction of cycle: #{direction_straight ? "#{index0[0]+1} -> #{index0[1]+1}" : "#{index0[1]+1} -> #{index0[0]+1}" }"
    i = index0[0]
    step = nil
    until i==index0[1]
      step = nil
      if direction_straight
	if basis[i][parents[i]]==1
	  step = [ x[i][parents[i]] - borders[i][parents[i]][0], [i,parents[i]] ] #обратная дуга
	  puts "step#{step[1][0]+1}#{step[1][1]+1} = x#{step[1][0]+1}#{step[1][1]+1} - d_low#{step[1][0]+1}#{step[1][1]+1}= #{x[step[1][0]][step[1][1]]} - #{borders[step[1][0]][step[1][1]][0]} = #{step[0]}"
	  step0 = step if step0[0] > step[0]
	elsif basis[parents[i]][i]==1
	  step = [ borders[parents[i]][i][1] - x[parents[i]][i], [parents[i],i] ]
	  puts "step#{step[1][0]+1}#{step[1][1]+1} = d_high#{step[1][0]+1}#{step[1][1]+1} - x#{step[1][0]+1}#{step[1][1]+1} = #{borders[step[1][0]][step[1][1]][1]} - #{x[step[1][0]][step[1][1]]} = #{step[0]}"
	  step0 = step if step0[0] > step[0]
	end
      else
	if basis[i][parents[i]]==1
	  step = [ borders[i][parents[i]][1] - x[i][parents[i]] , [i,parents[i]] ]
	  puts "step#{step[1][0]+1}#{step[1][1]+1} = d_high#{step[1][0]+1}#{step[1][1]+1} - x#{step[1][0]+1}#{step[1][1]+1} = #{borders[step[1][0]][step[1][1]][1]} - #{x[step[1][0]][step[1][1]]} = #{step[0]}"
	  step0 = step if step0[0] > step[0]
	elsif basis[parents[i]][i]==1
	  step = [ x[parents[i]][i] - borders[parents[i]][i][0] , [parents[i],i] ]
	  puts "step#{step[1][0]+1}#{step[1][1]+1} = x#{step[1][0]+1}#{step[1][1]+1} - d_low#{step[1][0]+1}#{step[1][1]+1}= #{x[step[1][0]][step[1][1]]} - #{borders[step[1][0]][step[1][1]][0]} = #{step[0]}"
	  step0 = step if step0[0] > step[0]
	end
      end
      #puts "step#{step[1][0]+1}#{step[1][1]+1} = #{step[0]}" unless step.nil?                    
      i = parents[i]
    end
    step = nil
    if direction_straight
	step = [ borders[index0[0]][index0[1]][1] - x[index0[0]][index0[1]], [index0[0],index0[1]] ]
	puts "step#{step[1][0]+1}#{step[1][1]+1} = d_high#{step[1][0]+1}#{step[1][1]+1} - x#{step[1][0]+1}#{step[1][1]+1} = #{borders[step[1][0]][step[1][1]][1]} - #{x[step[1][0]][step[1][1]]} = #{step[0]}"
	step0 = step if step0[0] > step[0]
    else
	step = [ x[index0[0]][index0[1]] - borders[index0[0]][index0[1]][0] , [index0[0],index0[1]] ]
	puts "step#{step[1][0]+1}#{step[1][1]+1} = x#{step[1][0]+1}#{step[1][1]+1} - d_low#{step[1][0]+1}#{step[1][1]+1}= #{x[step[1][0]][step[1][1]]} - #{borders[step[1][0]][step[1][1]][0]} = #{step[0]}"
	step0 = step if step0[0] > step[0]
    end
    #puts "step#{step[1][0]+1}#{step[1][1]+1} = #{step[0]}" unless step.nil?   
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
  
  def optimal? estimations, x, borders
    estimations.each_with_index.none? do |line, i| 
      line.each_with_index.any? do |estimate, j| 
	puts "Not optimal delta#{i+1}#{j+1} = #{estimate}" if !estimate.nil? && (x[i][j] == borders[i][j][0] ? estimate < 0: estimate > 0)
	estimate.nil? ? false : x[i][j] == borders[i][j][0] ? estimate < 0: estimate > 0
      end
    end
  end
  
  def choose_not_optimal_indexes basis, estimations, x, borders
    puts "Choose i0, j0 with the abs highest delta (not optimal)"
    min = []
    estimations.each_with_index do |line, i| 
      line.each_with_index do |edge, j| 
	min = [edge,[i,j]]if basis[i][j]==0 && !edge.nil? && (min.empty? || edge.abs > min[0].abs ) && (x[i][j] == borders[i][j][0] ? edge < 0: edge > 0)
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

  def find_init_basis trade, network, borders
    puts "First stage: "
    x_init, network_init, borders_init, cost_init, basis_init = *make_initial_parametres(trade, network, borders)
    print "Initial x: ", x_init
    print "Initial network: ", network_init, false
    print "Initial borders: ", borders_init, true
    #cost_init = make_cost_init trade
    print "Initial cost: ", cost_init
    #basis_init = deep_copy cost_init
    print "Initial basis: ", basis_init, false
    x = deep_copy x_init
    basis = deep_copy basis_init
    cost = deep_copy cost_init
    optimal = false
    counter = 1
    until optimal
      puts "Iteration ##{counter}"
      counter += 1
      x, basis, optimal = *iterate(basis, cost, network_init, x, borders_init)
      #raise "stop" if counter==8
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
  
  def make_initial_parametres trade, network, borders
    cost_init = Array.new(trade.size+1){[0]*(trade.size+1)}
    basis_init = Array.new(trade.size+1){[0]*(trade.size+1)}
    x_init = borders.map{|line| line.map{|item| item != 0  ? item[0] : 0 }}
    network_init = make_network_init trade, network
    borders_init = deep_copy borders
    closures = count_closures trade, x_init
    puts "Closures: "
    puts closures.each_with_index.map{|cl, i| "w#{i+1} = #{cl}"}
    x_init.each{|line| line << 0} << [0]*(x_init.size+1)
    borders_init.each{|line| line << 0} << [0]*(borders_init.size+1)
    closures.each_with_index do |tr, i|
      if tr >= 0 
	x_init[i][closures.size] = closures[i].abs
	network_init[i][closures.size] = 1
	borders_init[i][closures.size] = [0, closures[i].abs]
	cost_init[i][closures.size] = 1
	basis_init[i][closures.size] = 1
      else
	x_init[closures.size][i] = closures[i].abs
	network_init[closures.size][i] = 1
	borders_init[closures.size][i] = [0, closures[i].abs]
	cost_init[closures.size][i] = 1
	basis_init[closures.size][i] = 1
      end
    end
    [x_init, network_init, borders_init, cost_init, basis_init]
  end
  
  
  def count_closures trade, x_init
    trade.each_with_index.map{|tr, i| tr - x_init[i].inject{|s,x|s+x} + x_init.transpose[i].inject{|s,x|s+x}}
  end
  
  def make_cost_init trade
    trade.map do |cost|
      [0]*trade.size << (cost>= 0 ? 1 : 0)
    end << (trade.map{|cost| cost < 0 ? 1 : 0}+[0])
  end
  
  def make_network_init trade, network
    network.each_with_index.map do |point, i|
      point << 0
    end << ([0]*(trade.size+1))
  end
  
end
  
include PotentialMethodHelper

network = [
  [0,1,0,0,1,0,0,0],
  [0,0,1,0,0,0,1,0],
  [0,0,0,1,0,0,0,0],
  [0,0,0,0,0,1,0,0],
  [0,0,0,0,0,1,0,1],
  [0,0,0,0,0,0,0,1],
  [1,0,1,0,1,0,0,0],
  [0,0,1,1,0,0,1,0]
  ]
cost = [
  [0,5,0,0,9,0,0,0],
  [0,0,8,0,0,0,6,0],
  [0,0,0,13,0,0,0,0],
  [0,0,0,0,0,15,0,0],
  [0,0,0,0,0,4,0,1],
  [0,0,0,0,0,0,0,4],
  [9,0,21,0,25,0,0,0],
  [0,0,13,14,0,0,15,0]
  ]
borders = [
  [0,[15,50],0,0,[20,45],0,0,0],
  [0,0,[18,50],0,0,0,[20,55],0],
  [0,0,0,[10,47],0,0,0,0],
  [0,0,0,0,0,[4,23],0,0],
  [0,0,0,0,0,[5,24],0,[9,21]],
  [0,0,0,0,0,0,0,[12,37]],
  [[5,11],0,[8,46],0,[9,21],0,0,0],
  [0,0,[6,18],[8,17],0,0,[10,12],0]
  ]

trade = [70,46,-95,-35,-10,24,0,0]
	                        
x_init = [
  [0,50,0,0,55,0,0,0],
  [0,0,40,0,0,0,50,0],
  [0,0,0,37,0,0,0,0],
  [0,0,0,0,0,7,0,0],
  [0,0,0,0,0,13,0,15],
  [0,0,0,0,0,0,0,40],
  [10,0,32,0,18,0,0,0],
  [0,0,15,30,0,0,10,0]
  ]
	                        
basis_init = [
  [0,1,0,0,1,0,0,0],
  [0,0,1,0,0,0,1,0],
  [0,0,0,0,0,0,0,0],
  [0,0,0,0,0,0,0,0],
  [0,0,0,0,0,1,0,0],
  [0,0,0,0,0,0,0,0],
  [0,0,1,0,1,0,0,0],
  [0,0,0,0,0,0,0,0]
  ]

find_path trade, network, cost, borders#, x_init, basis_init
