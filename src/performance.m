% performance.m: This script computes precision and recall for any class 
% and show the confusion matrix.
%
% Input values:
%   Qns - N x S x K matrix (output of infer.m)
%   labels - N x K matrix (groud truth)
%
% added by
% Emanuele Sansone GCNU 15/12/14

function performance(Qns, labels)

[n K] = size(labels);
Qns = squeeze(sum(Qns,2));
Qns(find(Qns >= .5)) = 1;
Qns(find(Qns < .5)) = 0;
Qns = logical(Qns);
confusion = zeros(K,K);
for i = 1:K
   for j = 1:K
      if (i ~= j)
        confusion(i,j) = sum(Qns(:,i).*labels(:,j));  
      else
        confusion(i,i) = sum(Qns(:,i).*labels(:,i));
      end
   end
end
plotconfusion(labels',Qns')
accuracy = 0;
for i = 1:K
   precision = confusion(i,i)/sum(confusion(i,:))*100;
   recall = confusion(i,i)/sum(confusion(:,i))*100;
   fprintf('\nPrecision for class %d: %.2f%%', i, precision);
   fprintf('\nRecall for class %d: %.2f%%\n', i, recall);
   accuracy = accuracy + confusion(i,i);
end
accuracy = accuracy/n * 100;
fprintf('\nAccuracy: %.2f%%\n\n', accuracy);

