function[objf,iterations,maxPmis]=case_wo_thlim(str)
mpc = loadcase(str);
for p=1:length(mpc.branch(:,1))
    if mpc.branch(p,6)~=0||mpc.branch(p,7)~=0||mpc.branch(p,8)~=0
        mpc.branch(p,6)=0;
        mpc.branch(p,7)=0;
        mpc.branch(p,8)=0;
        
    end
end
overrides = struct('opf', struct('return_raw_der', 'g'));
mpopt = mpoption(overrides);
r=runopf(mpc,mpopt);
for p=1:length(r.bus(:,1))
    Pmis(p) = r.raw.g(p);
end
maxPmis = max(Pmis);

objf=r.f;
iterations = r.raw.output.iterations;
y=[objf,iterations];
fid = fopen(strcat('sol',str,'.csv'),'w');
fprintf(fid,'%6.5f\n',objf);
fprintf(fid,'%d\n',iterations);
fprintf(fid,'%12.8f',maxPmis);
end