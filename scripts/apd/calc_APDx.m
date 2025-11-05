%
%CALC APDx: calculates the action potential duration (APD) at x% of the
%repolarization phase
%
%Typical call(s): [vMIN, vMAX, MUV, APD90, vAPD90, tAPD90] = calc_APDx(dt, v, x)
%
%Input:
%  dt ... time discretization
%  v  ... vector containing the transmembrane potential (one action potential)
%  x  ... the desired value for repolarization, for instance, for APD50%
%  one should provide 50
%
%Output:
%  vMIN  ... minimum value of the transmenbrane potential
%  vMAX  ... maximum value of the transmenbrane potential
%  MUV   ... the maximum derivative of the transmenbrane potential
%  APDx  ... the action potential duration (APD) at x% of the repolarization
%  vAPDx ... transmembrane potential at x% of the repolarization
%  tAPDx ... time in which Vm = VAPDx
%
%Example: [vMIN, vMAX, MUV, APD90, vAPD90, tAPD90] = calc_APDx(dt, v, 90)
%
%(c) Fernando Campos 05-2008

function [vMIN, vMAX, MUV, APDx, vAPDx, tAPDx] = calc_APDx(dt, vm, x)

   % Time vector
   ti = 0;
   tf = (length(vm) - 1)*dt;
   t  = ti:dt:tf;

   % Time of the maximum upstroke velocity
   dvm         = gradient(vm, dt);
   [MUV, iMUV] = max(dvm);

   % Amplitude of the AP
   [vMAX, ivMAX] = max(vm);
   [vMIN, ivMIN] = min(vm);
   
   A = vMAX - vMIN;
   
   % The value of vm at x% of the repolarization
   vmx = (1 - x/100)*A;
   
   % Look for vmx in vm
   n     = length(vm);
   iAPDx = ivMAX;
   
   for i = ivMAX:1:n
      if (vm(i) <= min(vm) + vmx)
         iAPDx = i;
         break;
      end
   end
   % Calculating the APD at x% of repolarization
   APDx  = t(iAPDx) - t(iMUV);
   vAPDx = vm(iAPDx);
   tAPDx = t(iAPDx);
   
   disp(sprintf('\n'))
   disp(sprintf('min   = %.2f mV\n' , vMIN))
   disp(sprintf('max   = %.2f mV\n' , vMAX))
   disp(sprintf('APA   = %.2f mV\n' , vMAX - vMIN))
   disp(sprintf('MUV   = %.2f V/s\n', MUV))
   disp(sprintf('APD%d = %.2f ms\n' , x, APDx))
   disp(sprintf('LAT   = %.2f ms\n' , iMUV*dt))
end
