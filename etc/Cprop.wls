#!/usr/bin/env wolframscript -print all
If[Length[$ScriptCommandLine]<3,Quit[]];
p=ToExpression[$ScriptCommandLine[[2]]];
n=ToExpression[$ScriptCommandLine[[3]]];
p1=-p-1; (* change p1 to -p if rounding is not to the nearest *)
p2=2^p1;
d1m[e_]:=(1-e)/(1+e); (* \delta_1^- *)
d1p[e_]:=(1+e)/(1-e); (* \delta_1^+ *)
d2m[e_]:=((1-e)^(5/2))*Sqrt[1-(e*(2-e))/2]; (* \delta_2^- *)
d2p[e_]:=((1+e)^(5/2))*Sqrt[1+(e*(2+e))/2]; (* \delta_2^+ *)
dam[e_]:=(1-e)/d2p[e];
dap[e_]:=(1+e)/d2m[e];
fam[e_]:=(1-dam[e])/e;
fap[e_]:=(dap[e]-1)/e;
"fam["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[fam[p2]],n]]
"fap["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[fap[p2]],n]]
ddm[e_]:=d2m[e]/(1+e);
ddp[e_]:=d2p[e]/(1-e);
edm[e_]:=(ddm[e])^2-1;
edp[e_]:=(ddp[e])^2-1;
dfm[e_]:=(ddm[e]*(1-e))/(Sqrt[1+edp[e]]*((1+e)^(5/2)));
dfp[e_]:=(ddp[e]*(1+e))/(Sqrt[1+edm[e]]*((1-e)^(5/2)));
fem[e_]:=(1-dfm[e])/e;
fep[e_]:=(dfp[e]-1)/e;
"fem["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[fem[p2]],n]]
"fep["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[fep[p2]],n]]
efm[e_]:=dfm[e]^2-1;
efp[e_]:=dfp[e]^2-1;
dcm[e_]:=d1m[e]/(Sqrt[1+efp[e]]*Sqrt[1+e]);
dcp[e_]:=d1p[e]/(Sqrt[1+efm[e]]*Sqrt[1-e]);
cem[e_]:=(1-dcm[e])/e;
cep[e_]:=(dcp[e]-1)/e;
"cem["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[cem[p2]],n]]
"cep["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[cep[p2]],n]]
pel[e_]:=((1+2*d2p[e])*(1+e)+1)*(1+e);
"pel["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[pel[p2]],n]]
eps[e_]:=FullSimplify[e*(2+dap[e])+(e^2)*(1+dap[e])];
"eps["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[eps[p2]/p2],n]]
epp[e_]:=FullSimplify[Sqrt[2]*(eps[e]*(1+e)+e)];
"epp["<>ToString[p1]<>","<>ToString[n]<>"]="<>ToString[N[FullSimplify[epp[p2]/p2],n]]
hes[e_]:=FullSimplify[(dap[e]+(1+dap[e])*(1+e))*(1+e)];
(* TextForm[hes[e]] *)
he4[e_]:=FullSimplify[4-hes[e]];
TextForm[NSolve[he4[e]==0,e,Reals,WorkingPrecision->n]]
