function data = T6_full_data(withRadicals)
%T6_FULL_DATA Exact algebraic T6: 
%   Requires Symbolic Math Toolbox. No other files are needed.
%   DATA = T6_FULL_DATA() returns exact radicals and the six matrices.
%   DATA = T6_FULL_DATA(false) omits algebraic radical evaluation.
%
%       structure layout:
%       data.P             3-by-2 symbolic base matrix
%       data.C{i}          3-by-2 symbolic increment C_i (cell array)
%       data.prevertices{i} = data.P + C_1 + ... + C_{i-1}
%       data.kappa(i)      exact LARGER quadratic root, greater than 1
%       data.X{i}          data.prevertices{i} + data.kappa(i)*data.C{i}
%       data.states        rows [v_i, u_i, h_i, H_i], with h_i = -p(v_i)
%
%   Every rational coefficient is entered as a quoted sym string.
%   Do not round kappa or convert the matrices to double for exact checks.
%   The C and prevertices cells are also available withRadicals = false.
%   Each quadratic here has exactly one positive root (the larger root).

if nargin < 1, withRadicals = true; end
if ~islogical(withRadicals) || ~isscalar(withRadicals)
    error('T6:Input','withRadicals must be a logical scalar.');
end

data.s = [sym('765199014123/500000000000'); ...
    sym('281160055171/500000000000'); ...
    sym('563320140417/1000000000000'); ...
    sym('189450484127/125000000000'); ...
    sym('577607220591/1000000000000'); ...
    sym('1666793914751/1000000000000')];
data.g = [sym('-2602836919/100000000000'); ...
    sym('-624510145841092415937488516219915141/526438472020939525414833000000000000'); ...
    sym('546874620027/1000000000000'); ...
    sym('-815614411547431773606778729148371/14412701186397989519939000000000000'); ...
    sym('660009417549/1000000000000'); ...
    sym('103576165155174317635325259688127/1669854405838109370096150000000000')];
data.d = [sym('-6630126226357506810666297/100000000000000000000000000'); ...
    sym('353817231984652749549503401406967946618837466638109/526438472020939525414833000000000000000000000000000'); ...
    sym('-684975008733137960046293921/2000000000000000000000000000'); ...
    sym('-422432094116607657697224764155329005845095264777/3603175296599497379984750000000000000000000000000'); ...
    sym('-619945082853452890721477721/2000000000000000000000000000'); ...
    sym('547391871257466991839686648646005291057987633837/3339708811676218740192300000000000000000000000000')];
data.P = [sym('8061428053957/1000000000000'), sym('-197899143917/200000000000'); ...
    sym('-197899143917/200000000000'), sym('586061930617/500000000000'); ...
    sym('222840999665977/400000000000000'), sym('-1159812183363277/1000000000000000')];

% Primitive integer coefficients [A,B,C], with A > 0.
data.polyABC = [sym('3035408331540464108699611003589329043631777712566648387'), sym('-6929929572582847928641037403083595755327316750000000000'), sym('-799810989354915137500000000000000000000000000000000'); ...
    sym('8668418219873228467592807878623809567698096762435086746111332125639768295902115974626450678521635274489891'), sym('-51110774056796778537785917456980124084186943953254145583880650689543693621160863941664394679152444008264860'), sym('-107953940103462583416301269441326663900885208242641068552356993338708372376779962720967083533304020495700'); ...
    sym('14816192558498220545363847643078120345867515958884271591858897929342280228966248685530388742131692933621753'), sym('129900450549364132988835559679130772750633773252491593855075264329190330954909660485363693133242090328810374'), sym('-340402478216216101228875286782381129443117859192493599250573404457941069581084223438071288673312902034165352'); ...
    sym('2825125472179693735118045353932034648570666752397745650811051882487600938729387749165440524828375446581990871667400738712751616'), sym('5801499918333865802085209404703438532226998049010463141746172806194910033382694637069813764023870603240414534060603810878275320'), sym('-178922828532196595750363852958113185098127250912857346758735713447015242430740274238440540100289168911342290272395649314591089025'); ...
    sym('93630061952423393552639350962345697789661961494847403629663914991796674644404275460819829694512383454559'), sym('655189716571211057600650293288907934234023757114677819393460717069341260257597195653103461264925900518122'), sym('-749584443912103319097788095445573196007473302061676067400328324678063937500146311677282057557447348421081'); ...
    sym('49678145334751925171049243410591217436856127004047912145041238803382670236267961582029497754977338879'), sym('-47768266026739338045512434311384818934533398262274110815131576261542663740905860173632581259954677758'), sym('-1911663471172169861246127985798909959468958630377110943009231542315006495362101408396916495022661121')];

% Exact rational isolating bounds for the LARGER root of each quadratic.
data.kappaBounds = [sym('114157291648753693069854719964636963616320721176621962167794676627/50000000000000000000000000000000000000000000000000000000000000000'), sym('45662916659501477227941887985854785446528288470648784867117870651/20000000000000000000000000000000000000000000000000000000000000000'); ...
    sym('589831677008160575179206490341014846061213595474779686310991279867/100000000000000000000000000000000000000000000000000000000000000000'), sym('147457919252040143794801622585253711515303398868694921577747819967/25000000000000000000000000000000000000000000000000000000000000000'); ...
    sym('105590764775173524751613857520256441870480461418241211062144781179/50000000000000000000000000000000000000000000000000000000000000000'), sym('211181529550347049503227715040512883740960922836482422124289562359/100000000000000000000000000000000000000000000000000000000000000000'); ...
    sym('699737971561377649380139193896081287158609280289258413103152229007/100000000000000000000000000000000000000000000000000000000000000000'), sym('43733623222586103086258699618505080447413080018078650818947014313/6250000000000000000000000000000000000000000000000000000000000000'); ...
    sym('50045378859903078282745421789674882285825544245675956252623919201/50000000000000000000000000000000000000000000000000000000000000000'), sym('100090757719806156565490843579349764571651088491351912505247838403/100000000000000000000000000000000000000000000000000000000000000000'); ...
    sym('100003458368013388394934523710772641939646271801724887104535320791/100000000000000000000000000000000000000000000000000000000000000000'), sym('12500432296001673549366815463846580242455783975215610888066915099/12500000000000000000000000000000000000000000000000000000000000000')];

% Additional data for this configuration, not the old sloping barrier.
% The horizontal line u = -99/100 lies below every state and above the
% second Hugoniot branch from X_4 at v = v_5.
data.separating_line = sym('-99/100');
data.hprime = [sym(1);sym(100000);sym(1);sym('1/1000000');sym(100000);sym(1)];

% Symmetry-space coordinates are [v; u; E; h; Q].
% The corresponding 3-by-2 matrix is [v,u; u,h; E,Q].
data.pS = [data.P(1,1);data.P(2,1);data.P(3,1);data.P(2,2);data.P(3,2)];
data.CS = sym(zeros(5,6)); data.PS = sym(zeros(5,6));
data.C = cell(6,1); data.prevertices = cell(6,1);
prefix = data.pS;
for i=1:6
    g=data.g(i); d=data.d(i); s=data.s(i);
    data.CS(:,i)=[g;-g*s;d;g*s^2;-d*s];
    data.PS(:,i)=prefix; prefix=prefix+data.CS(:,i);
    z=data.CS(:,i); data.C{i}=[z(1),z(2);z(2),z(4);z(3),z(5)];
    z=data.PS(:,i); data.prevertices{i}=[z(1),z(2);z(2),z(4);z(3),z(5)];
end

% Check closure on the actual cell matrices, using only exact arithmetic.
data.closure = sym(zeros(3,2));
for i=1:6
    data.closure = data.closure + data.C{i};
end
data.closure = simplify(data.closure);
assert(isequal(data.closure,sym(zeros(3,2))),'Exact closure failed.');
assert(isequal(simplify(prefix-data.pS),sym(zeros(5,1))), ...
    'Exact prevertex closure failed.');

% Rational interval checks certify the larger root and kappa_i > 1.
for i=1:6
    a=data.polyABC(i,1); b=data.polyABC(i,2); c=data.polyABC(i,3);
    lo=data.kappaBounds(i,1); hi=data.kappaBounds(i,2);
    % Check the stored polynomial against the ACTUAL increments/prevertex.
    gi=data.g(i); si=data.s(i); di=data.d(i);
    u0=data.PS(2,i); h0=data.PS(4,i); Q0=data.PS(5,i);
    aa=gi^2*si^3;
    bb=-di*si-u0*gi*si^2+gi*h0*si;
    cc=Q0-u0*h0;
    assert(isAlways(gi~=0 & aa*b==bb*a & aa*c==cc*a,'Unknown','error'), ...
        'The stored quadratic does not match the actual flux equation.');
    assert(isAlways(a>0 & b^2-4*a*c>0,'Unknown','error'), ...
        'Invalid quadratic.');
    assert(isAlways(lo>1 & hi>lo & 2*a*lo+b>0,'Unknown','error'), ...
        'Invalid larger-root interval.');
    assert(isAlways(a*lo^2+b*lo+c<0 & a*hi^2+b*hi+c>0,'Unknown','error'), ...
        'Root interval does not bracket the larger root.');
end


% Certified velocity ordering and shock separation, using only rationals.
% Bounds use the coordinate order [v;u;E;h;Q], just like data.XS.
data.XSLower=sym(zeros(5,6)); data.XSUpper=sym(zeros(5,6));
for i=1:6
    for row=1:5
        c=data.CS(row,i); p=data.PS(row,i);
        if isAlways(c>=0,'Unknown','error')
            data.XSLower(row,i)=p+data.kappaBounds(i,1)*c;
            data.XSUpper(row,i)=p+data.kappaBounds(i,2)*c;
        else
            data.XSLower(row,i)=p+data.kappaBounds(i,2)*c;
            data.XSUpper(row,i)=p+data.kappaBounds(i,1)*c;
        end
    end
end
others=[1,3,4,5,6];
data.u2GapLower=data.XSLower(2,2)-data.XSUpper(2,others);
assert(all(isAlways(data.u2GapLower>0,'Unknown','error')), ...
    'X_2 is not strictly higher than all five other states.');
data.stateLineGapLower=data.XSLower(2,:)-data.separating_line;
assert(all(isAlways(data.stateLineGapLower>0,'Unknown','error')), ...
    'A state is not above the separating line.');
Llower=data.XSLower(1,5)-data.XSUpper(1,4);
Dlower=data.XSLower(4,5)-data.XSUpper(4,4);
u4lower=data.XSLower(2,4)-data.separating_line;
u4upper=data.XSUpper(2,4)-data.separating_line;
assert(isAlways(Llower>0 & Dlower>0 & u4lower>0,'Unknown','error'));
data.shockSquaredGapLower=Llower*Dlower-u4upper^2;
assert(isAlways(data.shockSquaredGapLower>0,'Unknown','error'), ...
    'The second Hugoniot branch is not strictly below the separating line.');

if withRadicals
    data.kappa=sym(zeros(6,1)); data.XS=sym(zeros(5,6));
    data.X=cell(6,1);
    for i=1:6
        a=data.polyABC(i,1); b=data.polyABC(i,2); c=data.polyABC(i,3);
        data.kappa(i)=(-b+sqrt(b^2-4*a*c))/(2*a);
        data.XS(:,i)=data.PS(:,i)+data.kappa(i)*data.CS(:,i);
        z=data.XS(:,i); data.X{i}=[z(1),z(2);z(2),z(4);z(3),z(5)];
    end
    data.v=data.XS(1,:).'; data.u=data.XS(2,:).';
    data.E=data.XS(3,:).'; data.h=data.XS(4,:).';
    data.H=data.E-data.u.^2/2;
    data.states=[data.v,data.u,data.h,data.H];
end
end
