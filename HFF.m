function [outputArg1] = HFF(inputArg1,inputArg2)
%HFF Summary of this function goes here
%   Detailed explanation goes here
[height, length, colors]=size(inputArg1);

partition_height = height/inputArg2;
partition_length = length/inputArg2;

Partitions = zeros(inputArg2, inputArg2, colors);
expanded_partitions = zeros(height, length, colors);

%split image into partitions, three color dimensions
for c = 1:colors
   for i = 1:inputArg2
      for j = 1:inputArg2
          sum_of_partition = 0;
          for p_i = 1:partition_height
              for p_j = 1:partition_length
                  i_val = (i-1)*partition_height+p_i;
                  j_val = (j-1)*partition_length+p_j;
                  sum_of_partition = sum_of_partition + inputArg1(i_val, j_val, c);
              end
          end
          
          u = sum_of_partition/(partition_height*partition_length);
          
          varriance = 0;
          
          for p_i = 1:partition_height
              for p_j = 1:partition_length
                  
                  difference = inputArg1(i_val, j_val, c) - u;
                  varriance = varriance + difference * difference;
              end
          end
          
          Partitions(i, j, c) = varriance;
          
      end
   end
end

for c = 1:colors
   for i = 1:inputArg2
      for j = 1:inputArg2
          for p_i = 1:partition_height
              for p_j = 1:partition_length
                  i_val = (i-1)*partition_height+p_i;
                  j_val = (j-1)*partition_length+p_j;
                  expanded_partitions(i_val, j_val, c) = Partitions(i, j, c); 
              end
          end
      end
   end
end

outputArg1 = expanded_partitions;


%focus measures, check if 0

%choose max

%weighted average
% a1 = a/(a+b)
% b1 = b/(b+b)

%output probibility

end

