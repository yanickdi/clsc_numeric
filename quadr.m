(* ::Package:: *)

ClearAll["Global`*"]
Print["Setting up equations..."]
qn[pn_, pr_] := 1 - (pn - pr) / (1 - delta)
qr[pn_, pr_] := (pn - pr) / (1 - delta) - pr /delta

(* Setting pn and rho to the case 1 \[Rule] rho \[GreaterEqual] 1*)
pn[wn_, pr_] := 1/2 * (1-delta+pr+wn)
rho[wn_, pr_] := (tau^(1/3) (-1+delta-pr+wn)^(2/3))/(2 (a (1-delta))^(1/3))
(* Retailer decided, now define manufacturer profit *)
profitm[wn_, pr_]:= qn[pn[wn,pr], pr]*(wn*(1- tau/rho[wn,pr]) - cn) + qr[pn[wn,pr],pr]*(pr - cr) + ((tau/rho[wn,pr])*qn[pn[wn,pr],pr] - qr[pn[wn,pr],pr])*s
(* The manufacturer lagrange function in the general form *)
Lm[wn_,pr_,mu1_,mu2_]:=profitm[wn,pr]-mu1*(-qr[pn[wn, pr], pr])-mu2*(qr[pn[wn, pr], pr]-tau/rho[wn, pr]*qn[pn[wn, pr], pr])

(* case 1a, rho \[GreaterEqual] 1 and qr = 0 *)
equations1a = {D[Lm[wn, pr, mu1, 0], wn] == 0, D[Lm[wn, pr, mu1, 0], pr] == 0, qr[pn[wn,pr], pr] == 0}

(* case 1b, rho \[GreaterEqual] 1 and 0 \[GreaterEqual] qr \[GreaterEqual] tau/rho*qn *)
equations1b = {D[Lm[wn, pr, 0, 0], wn] == 0, D[Lm[wn, pr, 0, 0], pr] == 0}

(* case 1c, rho \[GreaterEqual] 1 and qr = tau/rho*qn *)
equations1c = {D[Lm[wn, pr, 0, mu2], wn] == 0, D[Lm[wn, pr, 0, mu2], pr] == 0, qr[pn[wn,pr], pr] == tau/rho[wn,pr] * s}

(* Now case 2 - rho equals 1 *)
(* Setting pn and rho to the case 1 \[Rule] rho \[GreaterEqual] 1*)
pn[wn_, pr_] := 1/2 * (1-delta+pr+wn)
rho[wn_, pr_] := 1
(* Retailer decided, now define manufacturer profit *)
profitm[wn_, pr_]:= qn[pn[wn,pr], pr]*(wn*(1- tau/rho[wn,pr]) - cn) + qr[pn[wn,pr],pr]*(pr - cr) + ((tau/rho[wn,pr])*qn[pn[wn,pr],pr] - qr[pn[wn,pr],pr])*s
(* The manufacturer lagrange function in the general form *)
Lm[wn_,pr_,mu1_,mu2_]:=profitm[wn,pr]-mu1*(-qr[pn[wn, pr], pr])-mu2*(qr[pn[wn, pr], pr]-tau/rho[wn, pr]*qn[pn[wn, pr], pr])

(* case 2a, rho = 1 and qr = 0 *)
equations2a = {D[Lm[wn, pr, mu1, 0], wn] == 0, D[Lm[wn, pr, mu1, 0], pr] == 0, qr[pn[wn,pr], pr] == 0}

(* case 2b, rho \[GreaterEqual] 1 and 0 \[GreaterEqual] qr \[GreaterEqual] tau/rho*qn *)
equations2b = {D[Lm[wn, pr, 0, 0], wn] == 0, D[Lm[wn, pr, 0, 0], pr] == 0}

(* case 2c, rho \[GreaterEqual] 1 and qr = tau/rho*qn *)
equations2c = {D[Lm[wn, pr, 0, mu2], wn] == 0, D[Lm[wn, pr, 0, mu2], pr] == 0, qr[pn[wn,pr], pr] == tau/rho[wn,pr] * s}
Print[qr[pn[wn,pr], pr] == tau/rho[wn,pr] * s]
Print[D[Lm[wn, pr, 0, mu2], mu2]==0]

Print["Done, starting calculations..."]

i = 0
Do[
  param = {tau->tauVal, a->aVal, s->sVal, cr->crVal, cn->cnVal, delta->deltaVal};
  Print["--"]
  Print[param];

  Print["Case 1a"];
  Print[NSolve[equations1a/.param, {wn, pr,mu1}, Reals]];
  Print["Case 1b"];
  Print[NSolve[equations1b/.param, {wn, pr}, Reals]];
  Print["Case 1c"];
  Print[NSolve[equations1c/.param, {wn, pr, mu2}, Reals]];
  Print["Case 2a"];
  Print[NSolve[equations2a/.param, {wn, pr,mu1}, Reals]];
  Print["Case 2b"];
  Print[NSolve[equations2b/.param, {wn, pr}, Reals]];
  Print["Case 2c"];
  Print[NSolve[equations2c/.param, {wn, pr, mu2}, Reals]];
  i++;
, {tauVal, 0, 1, 0.1}, {aVal, 0, 0.1, 0.01}, {sVal, 0, 1, 0.1}, {crVal, sVal, 1, 0.1}, {cnVal, crVal, 1, 0.1}, {deltaVal, crVal, 1, 0.1}]
Print[i]




*
