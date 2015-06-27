function [ x, y, w ] = rule39 ( )

%*****************************************************************************80
%
%% RULE39 returns the rule of degree 39.
%
%  Discussion:
%
%    Order 39 (282 pts)
%    1/6 data for 39-th order quadrature with 54 nodes.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 June 2014
%
%  Author:
%
%    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
%    This MATLAB version by John Burkardt.
%
%  Parameters:
%
%    Output, real X(*), Y(*), the coordinates of the nodes.
%
%    Output, real W(*), the weights.
%
  x = [ ...
     -0.51907510039262467181240723005295, ...
     -0.27085312542591789026981871557739, ...
     -0.62535084249163554801261088188681E-01, ...
     -0.44171982481354024257677278787610, ...
      0.00000000000000000000000000000000, ...
     -0.15015342269982087499746548552200, ...
      0.00000000000000000000000000000000, ...
     -0.73242413817247197416912819947887E-01, ...
     -0.61576635195673085768179561363787, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
     -0.48003862478268162304965734842567, ...
     -0.36508484113357361074797983209786, ...
      0.00000000000000000000000000000000, ...
     -0.56713187575703966770973366858629, ...
     -0.14031613163580384457468777152879, ...
     -0.25349195822799364362686756031925, ...
     -0.23248267893550702450813771192708, ...
     -0.52163635975054715145448041789043, ...
      0.00000000000000000000000000000000, ...
     -0.90690106293315290819021535664027, ...
     -0.12932342290720507011192674546810, ...
     -0.13784410422442617726158008200667, ...
     -0.66294103217865996758685557932486, ...
     -0.64122202646501331223433885902278, ...
     -0.81530008183876068918941172976357, ...
     -0.90238671492070698682084691210471, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
     -0.36521604921152885859291119415780, ...
     -0.81466971749275774599927221958440E-01, ...
     -0.50100475560027347803268966859046, ...
      0.00000000000000000000000000000000, ...
     -0.42538501367251799249997046103522, ...
     -0.40434782435372039643598964415201, ...
      0.00000000000000000000000000000000, ...
     -0.16442016253358942481700772352343, ...
     -0.27207964379670182797824599426077, ...
     -0.31259592531293370325932513075940, ...
      0.00000000000000000000000000000000, ...
     -0.23061575503872384396183209443032, ...
     -0.73174496781264133588838797582992, ...
     -0.55942199595539976634905924759128, ...
     -0.28187070033874590044167993901641, ...
     -0.63181925527793262592959928916708, ...
     -0.95845151882990183009015464856753, ...
     -0.83806929801903041855622721765482, ...
      0.00000000000000000000000000000000, ...
     -0.73014202636963239486058389001507, ...
     -0.38481648423499752372487549629974, ...
     -0.74786883044034466526785200148550, ...
      0.00000000000000000000000000000000, ...
     -0.82569754738333487928005814615627, ...
      0.00000000000000000000000000000000 ];
  y = [ ... ...
     -0.56977344334973953415009221378062, ...
     -0.11615238463710973900986144045138, ...
     -0.16841403059063582703799040319519, ...
     -0.57273896147926731246264811881990, ...
     -0.48177391174488216182747312872034, ...
     -0.21805432181592896862251556772349, ...
     -0.43015355484213642389086647596671, ...
     -0.55672301202474317753018700126839, ...
     -0.57361430082863865627675251023354, ...
     -0.71062437062144605538166336071450E-01, ...
     -0.57425388147699645590520477695232, ...
     -0.39157484051927385758800981076725, ...
     -0.38260125928736646468021583934304, ...
      0.76210815817498206981551181976317, ...
     -0.37974002964326184444303529558949, ...
     -0.47005865572418340661540287763308, ...
     -0.36610435880584296669017948490237, ...
     -0.51041551846753395228797796020113, ...
     -0.54750015447026086336476242420687, ...
      0.11000415018640989305561291676418E+01, ...
     -0.55212750394609514626346834644682, ...
     -0.31033539412300214897853090795097, ...
     -0.40217042416953081560570497781101, ...
     -0.45041551128090115382074323796237, ...
     -0.55393733579619755005613977870791, ...
     -0.50261462137952603540673793837763, ...
     -0.57241680160442647300899412288108, ...
      0.58655605653182625196365551350355, ...
     -0.26028645488038103161700300407634, ...
     -0.50737982602401753884163334078742, ...
     -0.52365228029452559568466482181594, ...
     -0.50930584971577307589941102218584, ...
     -0.35752411415512963812671499048932, ...
     -0.45339755648003318420528024155758, ...
     -0.30993099181146801236771669846609, ...
      0.10264156710616244871686286975762E+01, ...
     -0.57254647239454006820290004274071, ...
     -0.27970449294975765316042730395238, ...
     -0.57171692937562625460547538450204, ...
      0.42124024273454953036405961333276, ...
     -0.54949681464841745643044481184475, ...
     -0.50173031881501531356011446242568, ...
     -0.45905866008218877559480919998558, ...
     -0.44792709733285852771474443606843, ...
     -0.51472732382392613411706398642860, ...
     -0.57244321215608923333186877139555, ...
     -0.54626015892447713295423347103218, ...
      0.16065133195264123594739464348022, ...
     -0.57182062792860608241992084882251, ...
     -0.54913090515973984664975730120048, ...
     -0.54631988756965195283384974383402, ...
      0.87522865731065508967027909558189, ...
     -0.57117161215160068889104816203223, ...
      0.11420800543675486817879702813617E+01 ];
  w = [ ... ...
      0.10297002991753660501136103351341E-02, ...
      0.71409200454957424138919555702907E-02, ...
      0.84670476850581833366989493797193E-02, ...
      0.10061525828127805509276236358357E-02, ...
      0.25008518613742033230286581854602E-02, ...
      0.82379503874691016163177015741387E-02, ...
      0.31645234824978289216154005178556E-02, ...
      0.28454670499457543125480683962546E-02, ...
      0.90921496761479735112870080703174E-03, ...
      0.53162773859582871544697296720480E-02, ...
      0.55174868571100333634113037140252E-03, ...
      0.57658467066484687646578171852548E-02, ...
      0.64457871132908107545457606397354E-02, ...
      0.24406168935034987856598180879293E-02, ...
      0.58315103312031023476998033502050E-02, ...
      0.61853878639532303634906561313174E-02, ...
      0.74949719610539269221700456157038E-02, ...
      0.49000462324377305053951530096155E-02, ...
      0.29889444058930598950088528689227E-02, ...
      0.43996809067738625074189199968923E-03, ...
      0.12505146482203086295392343808630E-02, ...
      0.97825105755546869735655412908223E-02, ...
      0.81418004583935088155390428216798E-02, ...
      0.50503317417279426125202432815639E-02, ...
      0.24951091420867912577241746709087E-02, ...
      0.30478657588270378768815691124488E-02, ...
      0.64260046142092336339976252566979E-03, ...
      0.38189427195873096307259859473584E-02, ...
      0.50663628845451396008230346642253E-02, ...
      0.51764441731398308561678437452717E-02, ...
      0.49181464598226297253857270629541E-02, ...
      0.46465662734804209629119572618546E-02, ...
      0.44591898920056385785598751636219E-02, ...
      0.64136254490796000834798012697594E-02, ...
      0.79881644631115754408349534866859E-02, ...
      0.10476275740484627938272086913711E-02, ...
      0.14500024170436655618156678953819E-02, ...
      0.96194255594278304779191550064611E-02, ...
      0.15445701615508816200489059359044E-02, ...
      0.46888846948789891123342520895612E-02, ...
      0.37215235791213064892823990494283E-02, ...
      0.37005626494781645355545034117062E-02, ...
      0.57612742538281776709493292272587E-02, ...
      0.76066250692246841971988592588158E-02, ...
      0.43715876705429943949816550672880E-02, ...
      0.42564655054254956754287557583451E-03, ...
      0.20530208586099517375295895636913E-02, ...
      0.57438263789397392552513701846135E-02, ...
      0.11429571862723729149581721669044E-02, ...
      0.36214122455131093638488339463188E-02, ...
      0.26987232674585178812827363148936E-02, ...
      0.24397941805574716783727001505373E-02, ...
      0.10319996026161764514307891440984E-02, ...
      0.11509579298275706789740693353329E-03 ];

  return
end
