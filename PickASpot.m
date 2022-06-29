%PickASpot
function[answer] = PickASpot()

%Nspots = 14;

rng('shuffle');
%switch Nspots
   % case 14
        % return a letter A-M
        B = randi(14,[5,1]);
        C = {'A','B'...
            'C','D','E'...
            'F','G','H'...
            'I','J','K',...
            'L','M','N'};
        A = C(B);
  %  case 72
        % return a number 1-72
      %  A = randi(72,[5,1]);
  %  otherwise
     %   error('unknown Nspots, choose 16 or 72 to match with cardboard template');
%end

%rn = randperm(6);
rn = randi(2,[5,1]);
conc0 = {'1%','0.01%'};
Conc= conc0(rn);

answer=[A, Conc];
end

