#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

header = ARGF.gets
puts %w<bene score>.join("\t")

ARGF.each do |line|
  fields = line.chomp.split("\t")
  bene = fields.shift
  ccs = fields.select { |x| x == '3' }.length
  puts [bene, ccs].join("\t")
end
