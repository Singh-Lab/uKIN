# --------------------------------------------------------------------------- #
# uKIN - Using Knowledge In Networks
# --------------------------------------------------------------------------- #
#
# Author: Borislav Hristov (borislav@cs.princeton.edu)
#
# I. Input 
#
# There are three required input files:
# 1) a network file 
# 2) a prior knowledge file containing a list of nodes (genes) known to be
# disease associated, possibly with weights on them
# 3) a file of newly implicated genes, each with a weight
#
# !!! Plus a required path to a Matlab executable since uKIN needs Matlab to do
# some heavy lifting of big matrices. The path should be specified as:
# matlab=/path/to/matlab after the three input files.
# 
# Additionally, the user may provide:
# 5) a custom value of the restart parameter alpha controlling how much 
# influence each type of information (prior and new) has.
# Default is alpha = 0.5
# 6) a custom value of the diffusion parameter gamma controlling how far the prior
# knowledge spreads. Default is gamma = 1.
# 7) output prefix which is used in the beginning of the name of 
# the output file
#
# II. Output
#
# output_prefix_results.txt is written in the uKIN directory. The file
# contains a list of candidate genes ranked by how frequently they are
# visited as the guided walks reach the stationary distribution.
#
# III. How to run
#
# Simply issue:
# ruby uKIN.rb network_file.txt prior_knowledge.txt new_genes.txt matlab=/path/to/matlab {alpha=0.4 gamma=0.8 output_prefix=my_output}
#
# 
# Note that uKIN is implemented in Ruby but requires Matlab for the heavy matrix 
# operations. 
# You can simply install Ruby via: sudo apt-get install ruby-full
# Matlab is typically provided for free to academic institutions. You can also
# download a 30 day free trial from their website.
#
# IV. Input File Formats
#
# 1. Network file: each line specifies an edge, white space delimited:
# GENE_ID GENE_ID
#
# 2. Prior knowledge: list of genes
# GENE_ID
# GENE_ID
#
# 3. Newly implicated genes: each line is white space delimited
# GENE_ID WEIGHT
# GENE_ID WEIGHT
#
#
# --------------------------------------------------------------------------- #



# --------------------------------------------------------------------------- #
# Parse User Input
# --------------------------------------------------------------------------- #

def parse_user_input
  if ARGV.size < 4
    puts "Insufficient inputs provided."
    exit
  end
  
  $network       = ARGV[0]
  $prior_file    = ARGV[1]
  $new_file      = ARGV[2]
    
  unless ARGV[3][0..6] == "matlab="
    puts "need to provide a path to the matlab executable matlab=/path/to/my/matlab"
    exit
  end
  
  $matlab_exe    = ARGV[3][7..-1]

  $alpha         = 0.5
  $gamma         = 1.0
  $out_file      = "uKIN"
  
  # see if any additional user inputs are desired to overwrite the defaults
  (4..(ARGV.size - 1)).each do |i|
    index = ARGV[i].index("=")

    if index.nil?
      puts "Illegal input format. Please use param=value without white spaces if adding any non-required parameters, i.e alpha=0.4"
      exit
    end
    
    param = ARGV[i][0..(index - 1)]
    value = ARGV[i][(index+1)..-1]
    
    case param
    when "alpha"
      $alpha = value.to_f
    when "gamma"
      $gamma = value.to_f
    when "output_prefix"
      $out_file = value
    else
      puts "Wrong parameter specified. Please use one of the alpha=, gamma=, or output_prefix=."
      exit
    end
  end
end

def run_with_user_inputs
  read_network
  
  read_prior_knowledge
  
  read_new_information
  
  propagate_prior_knowledge
  
  perform_guided_random_walks
end



# --------------------------------------------------------------------------- #
# Read Input Files
# --------------------------------------------------------------------------- #

def read_network
  $g = RubGraph.new
  
  File.foreach($network) do |line|
    words = line.delete("\n").split()
    $g.add_edge(words[0], words[1])
  end
  
  # fix the order of nodes
  $nodes = $g.nodes.keys
  
  puts "Read network file: #{$network}"
  puts " number of nodes = #{$g.nodes.keys.size}"
  puts " number of edges = #{$g.num_edges}"
end


def read_prior_knowledge
  $prior_kn = {}
  
  File.foreach($prior_file) do |line|
    words = line.delete("\n").split()

    if words.size > 1
      if $nodes.include?(words[0])
        $prior_kn[words[0]] = words[1].to_f
      else
        puts "gene #{words[0]} is not in the network."
      end
      
    else
      if $nodes.include?(words[0])
        $prior_kn[words[0]] = 1
      else
        puts "gene #{words[0]} is not in the network."
      end
    end
  end
  
  puts "Read prior knowledge file: #{$prior_file}\n  with #{$prior_kn.size} genes present in the network.\n"
end


def read_new_information
  $new_info = {}
  
  File.foreach($new_file) do |line|
    words = line.delete("\n").split()
    
    if $nodes.include?(words[0])
      $new_info[words[0]] = words[1].to_f
    else
      puts "gene #{words[0]} is not in the network."
    end
  end
  
  puts "Read new information file: #{$new_file}\n  with #{$new_info.size} genes present in the network.\n"
end


# --------------------------------------------------------------------------- #
# HASH FUNCTIONS
# --------------------------------------------------------------------------- #

# Read Hash Generic: reads a text file containing two columns of key value pairs
def read_text_file_hash file_name, is_int = false
  hsh = {}
  File.foreach(file_name) do |line|
    words = line.delete("\n").split()
    hsh[words[0]] = words[1]
    hsh[words[0]] = words[1].to_i if is_int
  end
  hsh
end

# make a deep copy of a hash
def deep_copy_hash g
  h = {}
  g.each_pair do |k,v|
    v = Marshal.load(Marshal.dump(v))
    h[k] = v
  end
  h
end

# sort the hash
def sort_hash hsh
  arr = []
  hsh.values.uniq.sort.reverse.each do |v|
    hsh.map{ |k,v1| (v1==v) ? k : nil }.compact.each { |e| arr << e }
  end
  arr
end

# reads a pre-ordered hash of floating values
def read_preordered_hash_f file_name
  hsh, ordered_genes = {}, []
  File.foreach(file_name) do |line|
    words = line.delete("\n").split()
    hsh[words[0]] = words[1].to_f
    ordered_genes << words[0]
  end
  
  return [hsh, ordered_genes]
end

# normalizes a hash of values
def normalize_hash_vector hsh
  sum = 0
  hsh.each_key { |k| sum += hsh[k] }
  hsh.each_key { |k| hsh[k] = (hsh[k]/sum).round(12) }
  return hsh
end


# --------------------------------------------------------------------------- #
# Small helper functions
# --------------------------------------------------------------------------- #

# return the fraction of covered patients
def get_frac_covered patients, cover
  return 0 if patients.size == 0
  (patients.keys.count { |pat| !(patients[pat]&cover).empty? }.to_f/patients.size).round($d_prec)
end

# return the number of covered patients
def get_num_covered patients, cover
  patients.keys.count { |pat| !(patients[pat]&cover).empty? }
end

# return the number of newly covered patients by the specified gene
def newly_cov_pat gene, cov_pat
  new_cov = Set.new
  return new_cov unless $genes.has_key? gene
  $genes[gene].each { |pat| new_cov << pat unless cov_pat.include? pat }
  new_cov 
end

# Reads a Basic file and returns an array of its lines
def read_simple_file file
  arr = []
  File.foreach(file) do |line|
    arr << line.delete("\n")
  end
  return arr
end


# --------------------------------------------------------------------------- #
# Array Class
# --------------------------------------------------------------------------- #
class Array
  def sum
    self.inject(:+)
  end
  
  def avg
    (self.inject(:+).to_f/self.size).round($d_prec)
  end
  
  def std
    mu = self.avg
    n = self.size - (self.size > 25 ? 1 : 0)
    (Math.sqrt(self.map { |r| (r - mu) ** 2}.inject(:+).to_f/n)).round($d_prec)
  end
  
  def median
      sorted = self.sort
      len = sorted.length
      return (sorted[(len - 1) / 2] + sorted[len / 2]) / 2.0
  end
  
  def avg_column i
    self.map { |r| r[i]}.inject(:+).to_f/self.size
  end
  
  def std_column i
    mu = self.avg_column i
    n = self.size - (self.size > 25 ? 1 : 0)
    Math.sqrt(self.map { |r| (r[i]- mu) ** 2}.inject(:+).to_f/n)
  end

  def ratio_col i, j
    self.map { |r| r[i].to_f/r[j]}
  end
end


# --------------------------------------------------------------------------- #
# Simple Graph Class
# --------------------------------------------------------------------------- #
class RubGraph
  attr_reader :nodes
  
  def initialize(id = 'id_default')
    @id = id
    @nodes = {}
  end

  def add_node n
    @nodes[n] = {} unless @nodes.include? n
  end
  
  def add_edge e1, e2
    add_node e1; add_node e2
    @nodes[e1][e2] = 1 and @nodes[e2][e1] = 1
  end
  
  def remove_node n
    raise "node #{n} not in the graph" unless @nodes.has_key?(n)
    @nodes.delete n
    @nodes.each { |_, v| v.delete n }
  end
  
  def remove_edge n1, n2
    @nodes[n1].delete n2
    @nodes[n2].delete n1
  end
  
  def size
    @nodes.size
  end
  
  def degree v
    @nodes[v].size
  end
  
  def neighbors v
    @nodes[v].keys
  end
  
  def all_nodes_except v
    @nodes.keys.clone - ([]<<v)
  end
  
  def num_nodes
    @nodes.keys.size
  end
  
  def num_edges
    sum = 0
    @nodes.keys.each do |v|
      sum += @nodes[v].keys.size
    end
    sum/2
  end
end


# --------------------------------------------------------------------------- #
# The Pickup Library
# --------------------------------------------------------------------------- #

class Pickup
  attr_reader :list, :uniq
  attr_writer :pick_func, :key_func, :weight_func

  def initialize(list, opts={}, &block)
    @list = list
    @uniq = opts[:uniq] || false
    @pick_func = block if block_given?
    @key_func = opts[:key_func]
    @weight_func = opts[:weight_func]
  end

  def pick(count=1, opts={}, &block)
    func = block || pick_func
    key_func = opts[:key_func] || @key_func
    weight_func = opts[:weight_func] || @weight_func
    mlist = MappedList.new(list, func, uniq: uniq, key_func: key_func, weight_func: weight_func)
    result = mlist.random(count)
    count == 1 ? result.first : result
  end

  def pick_func
    @pick_func ||= begin
      Proc.new do |val|
        val
      end
    end
  end

  class CircleIterator
    attr_reader :func, :obj, :max, :key_func, :weight_func

    def initialize(obj, func, max, opts={})
      @obj = obj.dup
      @func = func
      @max = max
      @key_func = opts[:key_func] || key_func
      @weight_func = opts[:weight_func] || weight_func
    end

    def key_func
      @key_func ||= begin
        Proc.new do |item|
          item[0]
        end
      end
    end

    def weight_func
      @weight_func ||= begin
        Proc.new do |item|
          item[1]
        end
      end
    end

    def each
      until obj.empty?
        start = 0
        obj.each do |item|
          key = key_func.call(item)
          weight = weight_func.call(item)

          val = func.call(weight)
          start += val
          if yield([key, start, max])
            obj.delete key
            @max -= val
          end
        end
      end
    end
  end

  class MappedList
    attr_reader :list, :func, :uniq, :key_func, :weight_func

    def initialize(list, func, opts=nil)
      if Hash === opts
        @key_func = opts[:key_func]
        @weight_func = opts[:weight_func] || weight_func
        @uniq = opts[:uniq] || false
      else
        if !!opts == opts
          # If opts is explicitly provided as a boolean, show the deprecated warning.
          warn "[DEPRECATED] Passing uniq as a boolean to MappedList's initialize method is deprecated. Please use the opts hash instead."
        end

        @uniq = opts || false
      end

      @func = func
      @list = list
      @current_state = 0
    end

    def weight_func
      @weight_func ||= begin
        Proc.new do |item|
          item[1]
        end
      end
    end

    def each(&blk)
      CircleIterator.new(@list, func, max, key_func: @key_func, weight_func: weight_func).each do |item|
        if uniq
          true if yield item
        else
          nil while yield(item)
        end
      end
    end

    def random(count)
      raise "List is shorter then count of items you want to get" if uniq && list.size < count
      nums = count.times.map{ rand(max) }.sort
      return [] if max == 0
      get_random_items(nums)
    end

    def get_random_items(nums)
      current_num = nums.shift
      items = []
      each do |item, counter, mx|
        break unless current_num
        if counter%(mx+1) > current_num%mx
          items << item
          current_num = nums.shift
          true
        end
      end
      items
    end

    def max
      @max ||= begin
        max = 0
        list.each{ |item| max += func.call(weight_func.call(item)) }
        max
      end
    end
  end
end


# --------------------------------------------------------------------------- #
# uKIN Step1
# --------------------------------------------------------------------------- #

def propagate_prior_knowledge
  puts "First will spread the prior knowledge over the network..."
  f1 = File.open("temp/#{$out_file}_legend1.txt", 'w')
  f2 = File.open("temp/#{$out_file}_adj1.txt", 'w')
  f3 = File.open("temp/#{$out_file}_b1.txt", 'w')

  # 0. initialize the hash to transition to neighbors
  initialize_trans_prob_to_neighbors

  # 1. compute the value of each cell and write it down
  $nodes.each do |v1|
    f1.puts "#{v1}"  
    $nodes.each do |v2|
      v1_v2 = ($hsh_transit_prob[v2].has_key?(v1) ? $hsh_transit_prob[v2][v1] : 0)
      
      f2.print "#{(v1_v2).to_f.round(9)} "
    end
    f2.print "\n"
  end
  
  # 2. write down the querry vector b
  $nodes.each do |v1|
    f3.puts "#{($prior_kn.has_key?(v1) ? $prior_kn[v1] : 0)}"
  end
  f1.close; f2.close; f3.close
  
  # 3. write the matlab config to propagate the prior knowledge
  write_matlab_config_propagate_prior_knowledge
  
  # 4. propagate it
  execute_matlab_config("#{$out_file}_prior_config.m")
  
  # 5. normalize and sort the scores
  sort_and_normalzie_diffusion_scores
  
  puts "Done spreading prior knowledge. Initiating guided walks. This will take a while...\n"
end


def sort_and_normalzie_diffusion_scores
  hsh, i = {}, 0
  
  # 0. read the order of nodes
  arr_legend = read_simple_file("temp/#{$out_file}_legend1.txt")
  
  # 1. read the diffusion scores
  File.foreach("temp/#{$out_file}_ans1.txt") do |line|
    hsh[arr_legend[i]] = line.delete("\n").to_f + 0.00001
    i += 1
  end
  
  # 2. normalize them
  hsh = normalize_hash_vector(hsh)
  
  # 3. ready to use them
  $hsh_flow = hsh
end


def initialize_trans_prob_to_neighbors
  $hsh_transit_prob = Hash.new {|hash, key| hash[key] = {} }
  
  $g.nodes.each_key do |v|
    $g.neighbors(v).each { |e| $hsh_transit_prob[v][e] = 1 }
  end
end


# --------------------------------------------------------------------------- #
# uKIN Step2
# --------------------------------------------------------------------------- #

def perform_guided_random_walks
  # 0. make sure the starting frequencies sum to 1 and every node has epsilon
  normalize_starting_frequencies
  
  # 1. use the prior knowledge to bias the transition probabilities
  compute_transition_prob_from_prior_knowledge
  
  # 2. write down the transitional probability matrix P
  fo = File.open("temp/#{$out_file}_adj2.txt", 'w')
  
  $nodes.each do |v1|
    $nodes.each do |v2|
      v1_v2 = $hsh_transit_prob[v1].has_key?(v2) ? "#{$hsh_transit_prob[v1][v2]} " : "#{$alpha*$new_info[v2]} "
      fo.print "#{(v1_v2).to_f.round(10)} "
    end
    fo.print "\n"
  end
  fo.close
  
  # 3. write the matlab config to execute the guided walks
  write_matlab_config_guided_random_walk
  
  # 4. do it
  execute_matlab_config("#{$out_file}_guided_config.m")
  
  # 5. sort and write down the scores
  sort_write_scores
  
  # 6. do clean up of temporary files
  clean_up
  
  puts "Done performing guided walks and computing everything,\nwriting final scores to #{$out_file}_results.txt"
end


def normalize_starting_frequencies
  # give every node an epsilon 
  $nodes.each do |v|
    if $new_info.has_key?(v)
      $new_info[v] += 0.00001
    else
      $new_info[v] = 0.00001
    end
  end
  
  # and then normalize
  $new_info = normalize_hash_vector($new_info)
end


def compute_transition_prob_from_prior_knowledge
  # 1. define transitional probability matrix P
  $hsh_transit_prob = Hash.new {|hash, key| hash[key] = {} }

  # 2. for every node compute the transitional probabilities to its neighbors
  $g.nodes.each_key do |v|
    sum = 0
    $g.neighbors(v).each { |e| sum += $hsh_flow[e] }
    
    $g.neighbors(v).each { |e| $hsh_transit_prob[v][e] = $hsh_flow[e].to_f/sum }
  end
  
  # 3. Allow for restarts with coef alpha
  $g.nodes.each_key do |v|
    $g.neighbors(v).each { |e| $hsh_transit_prob[v][e] = ((1-$alpha)*$hsh_transit_prob[v][e] + $alpha*$new_info[e]) }
  end
end


def sort_write_scores
  $hsh_visit_prob, i = {}, 0
  
  # 0. read the order of nodes
  arr_legend = read_simple_file("temp/#{$out_file}_legend1.txt")
  
  # 1. read the eigevalues computed by matlab
  File.foreach("temp/#{$out_file}_ans2.txt") do |line|
    v = line.delete("\n").split()[0].to_f.abs
    $hsh_visit_prob[arr_legend[i]] = v
    i += 1
  end
  
  # 2. normalize the visit probabilities
  $hsh_visit_prob = normalize_hash_vector($hsh_visit_prob)
  
  # 3. sort the nodes based on the visit probability
  ordered_genes = sort_hash($hsh_visit_prob)
  
  # 4. write down the sorted visit probabilities
  fo = File.open("#{$out_file}_results.txt", 'w')
  ordered_genes.each { |e| fo.puts "#{e} #{$hsh_visit_prob[e].round(5)}" }
  fo.close
end


def clean_up
  FileUtils.rm("temp/#{$out_file}_legend1.txt")
  FileUtils.rm("temp/#{$out_file}_ans1.txt")
  FileUtils.rm("temp/#{$out_file}_ans2.txt")
  FileUtils.rm("temp/#{$out_file}_adj1.txt")
  FileUtils.rm("temp/#{$out_file}_adj2.txt")
  FileUtils.rm("temp/#{$out_file}_b1.txt")
  FileUtils.rm("temp/#{$out_file}_guided_config.m")
  FileUtils.rm("temp/#{$out_file}_prior_config.m")
end


# --------------------------------------------------------------------------- #
# Matlab
# --------------------------------------------------------------------------- #

def write_matlab_config_propagate_prior_knowledge
  fw = File.open("temp/#{$out_file}_prior_config.m", 'w')
  
  fw.puts "% read the input adjacency matrix"
  fw.puts "A = importdata(\"#{$matlab_dir}/#{$out_file}_adj1.txt\");"
  fw.puts "S = diag(sum(A));"
  fw.puts "I = eye(size(A,1));"
  fw.puts "y = #{$gamma};"
  fw.puts "b = importdata(\"#{$matlab_dir}/#{$out_file}_b1.txt\");"
  fw.puts ""
  fw.puts "L = -(A-S-y*I);"
  fw.puts "p = L^-1*b;"
  fw.puts ""
  fw.puts "% write the answer"
  fw.puts "dlmwrite(\"#{$matlab_dir}/#{$out_file}_ans1.txt\",p,\" \")"
  fw.close
end


def write_matlab_config_guided_random_walk
  fw = File.open("temp/#{$out_file}_guided_config.m", 'w')
  
  fw.puts "% read the transitional probabilities"
  fw.puts "mk = importdata(\"#{$matlab_dir}/#{$out_file}_adj2.txt\");"
  fw.puts "mk = transpose(mk);"
  fw.puts ""
  fw.puts "% compute the eigenvalues"
  fw.puts "[Vk,Dk] = eig(mk);"
  fw.puts ""
  fw.puts "% find where the 1 is"
  fw.puts "[row col] = ind2sub(size(Dk), find(Dk > 0.99999 & Dk < 1.00001))"
  fw.puts ""
  fw.puts "% in case of lower precision"
  fw.puts "if isempty(row)"
  fw.puts " [row col] = ind2sub(size(Dk), find(Dk > 0.99 & Dk < 1.01))"
  fw.puts "end"
  fw.puts ""
  fw.puts "% write the answer"
  fw.puts "dlmwrite(\"#{$matlab_dir}/#{$out_file}_ans2.txt\",Vk(:,col),\" \")"
  fw.close
end

# run matlab- expects a working version of Matlab given by $matlab_exe
def execute_matlab_config config_file
  cmd = "#{$matlab_exe}/matlab -nodisplay -nosplash -nodesktop -r \"run('temp/#{config_file}');exit;\""
  
  # run the command
  x = `#{cmd}`
end


# --------------------------------------------------------------------------- #
# Required libraries and basic paths and constants 
# --------------------------------------------------------------------------- #
require 'set'
require 'fileutils'

# Set the base directory to that containing this file
$BASE_DIR = (__dir__)
Dir.chdir $BASE_DIR
$matlab_dir = "#{$BASE_DIR}/temp"
$MIN_c, $d_prec = -0.000_001, 4



if __FILE__ == $PROGRAM_NAME
  
  parse_user_input
  run_with_user_inputs
end
