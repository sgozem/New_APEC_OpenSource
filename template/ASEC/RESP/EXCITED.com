%nproc=4
%mem=10000MB
#P cis=(Root=1,NStates=1)/6-31G* density=current SCF=Tight Pop=MK IOp(6/33=2) charge NoSymm

Gaussian template RESP charges

