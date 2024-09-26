function [LZ] = getMultiscaleLZ(M,scal)
%%%%transform an EEG signal into binary series at different scales
%to select scales: sampling rate/frec--------> alpha beta...

%%apply LZ to the binary series and store in a matrix


%M=EEG.data(1,:,1);


%%%scales. use odd number!!
%scal=[27 11 7];

for i=1:size(scal,2)
    if mod(scal(1,i),2)==0
        scal(1,i)=scal(1,i)+1;
    end

%%%side
sides=max(scal)-1;

%%%%%%%%%output
B=zeros(length(scal),length(M)-sides);
LZ=zeros(1,length(scal));


%%%%binary matrices
for j=1:length(scal)
     for i=sides/2:length(M)/2
         if M(i)> median(M(1+i-(scal(j)-1)/2:i+(scal(j)-1)/2))
        B(j,i)=1;
         end
     end
end

%%%LZ
for j=1:length(scal)
[C, H] = calc_lz_complexity(B(j,:), 'exhaustive', 1);
LZ(1,j)=C;
end

end

