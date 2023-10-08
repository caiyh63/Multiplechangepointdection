clear
fid  = fopen('ALTC.txt');
data = fscanf(fid, '%d');
fclose(fid) 

run ./MultipleChangePointDetection.m

TimeIndex=1979:2014

[ProbHypothesis, SampleOfHypothesisIndex, SampleOfChangePoint, SampleOfRate] = MultipleChangePointDetection(data, TimeIndex(1));


[AA,BB] = max(ProbHypothesis);

% 输出文件
fid1 = fopen('ProbHypothesis.txt','w');
fprintf(fid1, '%f\n', ProbHypothesis);
fclose(fid1) 
BB
if (BB > 1)
	zz = SampleOfChangePoint(1:BB-1,:,BB);
	zz = zz(:, zz(1,:)~=0);
end 
        
for p = 1:BB-1
	fid2 = fopen(['SampleOfChangePoint',num2str(p),'.txt'],'wt');
	fprintf(fid2, '%d\n',zz(p,:));
	fclose(fid2);
end  