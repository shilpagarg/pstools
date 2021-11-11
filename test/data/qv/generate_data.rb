fasta_str = String.new
1000.times.each do |n|
  genome_base = ['A', 'T', 'G', 'C'].sample
  fasta_str << genome_base
  # fasta_str << "\n" if (n + 1) % 60 == 0
end

puts fasta_str
