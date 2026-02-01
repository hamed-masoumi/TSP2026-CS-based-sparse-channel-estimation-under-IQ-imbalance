function  [DFTx, iDFTx] = func_ovcmp_dict(N,nx)


Gr = nx*N;
t1 = 0:(Gr-1);
t2 = 0:(N-1);
T=(t2.')*t1;

iDFTx = (1/sqrt(N))*exp(1i*2*pi*T/Gr); % "N x nxN" overcomplete idft-matrix-based dictionary

DFTx = (1/sqrt(N))*exp(-1i*2*pi*T/Gr); % "N x nxN" overcomplete idft-matrix-based dictionary

end