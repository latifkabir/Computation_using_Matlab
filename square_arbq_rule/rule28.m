function [ x, w ] = rule28 ( n )

%*****************************************************************************80
%
%% RULE28 returns the rule of degree 28.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 July 2014
%
%  Author:
%
%    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
%    This MATLAB version by John Burkardt.
%
%  Reference:
%
%    Hong Xiao, Zydrunas Gimbutas,
%    A numerical algorithm for the construction of efficient quadrature
%    rules in two and higher dimensions,
%    Computers and Mathematics with Applications,
%    Volume 59, 2010, pages 663-676.
%
%  Parameters:
%
%    Input, integer N, the number of nodes.
%
%    Output, real X(2,N), the coordinates of the nodes.
%
%    Output, real W(N), the weights.
%
  xs = [ ...
    0.2827621856762107,-.8043038962721508, ...
    -.3395780944817125,0.9664085101515743, ...
    -.1096329656209190E-01,0.8193786706976756, ...
    -.9286317198271657,0.4198022792026188, ...
    -.8604918125918429,0.8109783092388664, ...
    0.5614606479838853,0.1250139096754087E-01, ...
    0.6890874578865651,-.8003785261499256, ...
    -.7186772992712774,0.2097379198022608E-01, ...
    -.9904175004752394,0.5247683324321177, ...
    0.9879057318966782,0.9872995514825285, ...
    -.2010899760026548,-.7150754139309137, ...
    -.9744175518522036,-.2743957404149748, ...
    -.8801155850915744E-01,-.9120488211036881, ...
    0.5997911603089924E-01,-.9822903323396509, ...
    0.9933058426309856,-.1094425840746250, ...
    0.8016032845806299,-.3785896040093871, ...
    0.9908494000362051,-.9920099650417782, ...
    -.5334142633484058,-.6081971936243052, ...
    -.9913559247973139,-.4420258077139650, ...
    0.9904922021473850,-.2701059581013934, ...
    -.6762870118531359,0.4550203878192706, ...
    -.9933274734742550,0.1025031015215653, ...
    -.9908021888968861,0.4018136560549200E-02, ...
    0.1830465948969627,-.9190332253043175, ...
    -.7967293798867977,0.8910586530553062, ...
    0.9738677077447094,0.4922211332598602, ...
    -.9935974988322385,-.7469015595759197, ...
    -.7474851457426346,-.8461768034166064, ...
    0.7263311602661456,0.1838655693690960, ...
    -.8924925604346934,0.3473229569871696, ...
    0.7888084262939273,-.1622547625661196, ...
    0.1708950639909230,-.9858834062192522, ...
    0.8595474037303855,-.2487721441258377, ...
    0.9897369384324803,-.8592565754231390, ...
    0.3552890037180957,0.9343486754942414, ...
    -.5982760380368914,0.9012497934650560, ...
    -.8867898943082418,-.2973823530897202, ...
    0.2929780512625694,0.5109886162010564, ...
    -.6362138293509048,0.9840555620580882, ...
    0.6242255274464619,0.9500492760473710, ...
    0.6206677950945041,0.3667200892151265, ...
    -.3354257519635134,0.6532050712129780, ...
    0.8340712667685216,0.9517120425530485, ...
    0.4939764190001616,0.8519117390868273, ...
    -.4641713197529147,0.3096864797777032, ...
    -.9951406550164114,0.6849334791605266, ...
    0.5377487266764005,0.6596506837798173, ...
    -.6190982139272722,0.9676681958645598, ...
    -.8663734426507712,0.9229408505602535, ...
    0.8925298503917507,0.1376340175855351, ...
    -.4255923534036878,-.5069417238664446, ...
    -.9581264401591332,0.7956423904664977, ...
    -.9438198089736733,-.9602193183322451, ...
    -.9546324794056378,0.7865853669573317, ...
    0.7640178616435913,-.9037187541181603E-01, ...
    0.6583886775017005,0.9394184495408773, ...
    0.1210974568433837,0.3181846883529323, ...
    -.9476649688000099,-.7447346202738520, ...
    0.7669083321105266,0.3177180418892223, ...
    0.9005374034994568,0.4905235558364557, ...
    -.2799102574207112,-.7756668408683119, ...
    -.7847585523489610,-.5947448678197160, ...
    -.6835954896209653E-01,-.4603444798372620, ...
    0.1267789889405305,0.9957007235251634, ...
    -.9555617449567378,0.4926175951876377, ...
    0.9580744421864001,0.6544409855907934, ...
    -.2589485669037206,-.4625383187120728, ...
    0.8779137327105354,-.8753270142023138, ...
    -.7680237343653598,-.4646141113930034, ...
    -.6237643000374615,-.1004225708457810, ...
    0.2934408497776114,-.8838031019815079, ...
    -.8252198025376763E-01,0.1013709068378068, ...
    -.6355163316600663,-.2814473462899598, ...
    0.9991999747547105 ];
  ys = [ ...
    -.2146915742023843,-.3602259227928870, ...
    -.2366070789246893,-.1067737835329393, ...
    -.7832467667692294E-01,0.5090533388398083E-02, ...
    0.5655320173536176E-02,-.9091460408400718, ...
    -.3610436689009546,0.1173540113566495, ...
    -.7580822207698504,0.7636857906206649, ...
    -.6742700466367159,0.3605821467647248, ...
    0.2565264593758684,0.8682397249203467, ...
    -.9887795883901583,-.8777099834380453, ...
    -.9937283708213308,0.9906296802548810, ...
    0.8788092453036578,0.9706849956765652, ...
    -.7829696251737753E-01,0.9937649402523892, ...
    -.1524877190213388,-.2242135513466236, ...
    0.9919534212785273,0.9960588898544035, ...
    -.5034894383669408,0.9599951355029216, ...
    -.5645937597081716,0.9447996831444674, ...
    0.2711989824047218,0.1447359435070455, ...
    0.9871535260388329,0.9140753757822478, ...
    -.2960350108035981,0.8371402678277077, ...
    -.7650309750102933,-.3188681363434392, ...
    -.9938368658433656,-.6339545710762906, ...
    -.6212095493273496,-.9934382054756534, ...
    0.4715389357436094,-.7215460671816014, ...
    0.6773640866690520,0.9779955899159012, ...
    0.9970513707078691,-.1551145670991118, ...
    -.3024178806712221,-.9926567100374114, ...
    -.8660993570031557,-.4956056862948234, ...
    0.8382998942855415,0.9259687805728221, ...
    -.8351882927550272,0.9386358290968080, ...
    -.9900384261478909,-.4464745325070345, ...
    -.9936162118016632,-.8235872679034352, ...
    -.5983594352543780,0.9398154839397104, ...
    -.7371773517980565,0.7366500643247871, ...
    0.8982625204679442,0.1211489804261566, ...
    0.9859353129167986,-.8465334600907201, ...
    0.7342041846853576,-.4349009915613299, ...
    0.5008348384608579,-.9942925485328575, ...
    -.9626593056228334,-.2709700923096490, ...
    0.5944472663027801E-01,-.9248481403512147, ...
    -.4757878778805171,0.9921234070894396E-01, ...
    0.9864359054124004,0.5564466223050852, ...
    -.9075561035062022,0.6130603473347437, ...
    -.9196402734279348,-.6285797630023995, ...
    0.7360652264674404,0.9946281004658378, ...
    -.1056502846403645,0.8343858336015567, ...
    0.7795434558480316,0.2437361149690468, ...
    0.4083791315159580,-.8237952273107227E-01, ...
    -.6450626080823243,0.4764961946576270, ...
    0.7429354103100736,-.9752458082257867, ...
    0.3137605530793695,0.7602258278505522E-01, ...
    0.6028197206610887,-.9665994542595457, ...
    -.9448687562971092,0.7450892979512898, ...
    0.8580812769827701,0.6469685322672528, ...
    -.4730667721016463,0.4691067254414554, ...
    0.9468223400568799,0.5747271977509141, ...
    0.8487224420002684,0.9604420521431418, ...
    0.4325756220233750,-.1154193690840114, ...
    0.2996859346038543,0.6079262704780987, ...
    -.3014753479396839,0.2694764821728903, ...
    0.6210587156046703,0.9250623306784818, ...
    0.4091950530277537,-.1149400908487048, ...
    -.9509596164554358,0.4489210023756193, ...
    0.2491367754635779,0.2328542094943754, ...
    -.3147181280829211,0.6575404381717325, ...
    -.7572545223612748,0.8744761196977546E-01, ...
    0.7809511457782177,-.9573612252124380, ...
    0.5984333435621010E-01,-.4795365437914892, ...
    0.8758201084343584,-.6288667697841542, ...
    -.7704960778314179,-.7807547746986493, ...
    -.2948221971565224,-.9595773298541288, ...
    -.7859464731741874,-.8701241408788146, ...
    -.4848613986697829,-.8903571922698699, ...
    -.8832911479917025,-.6435327316304572, ...
    -.5075800418230667E-01 ];
  ws = [ ...
    0.1740144063407695E-02,0.6784563755290567E-02, ...
    0.1209466540188899E-01,0.5276845120506419E-02, ...
    0.2120637597472411E-01,0.1500950169216209E-01, ...
    0.1016605179698056E-01,0.1052237804962158E-01, ...
    0.1389531092331450E-01,0.1799231574573814E-01, ...
    0.1681749672744976E-01,0.2190516172058484E-01, ...
    0.1802050074136916E-01,0.1929827700816552E-01, ...
    0.2322655354121691E-01,0.1715984331158577E-01, ...
    0.7585152723116851E-03,0.1423586238543194E-01, ...
    0.6209204642340611E-03,0.8055582944619972E-03, ...
    0.1686431398199053E-01,0.6130199135277467E-02, ...
    0.8266675752855135E-02,0.3947441053087870E-02, ...
    0.3763432904885534E-01,0.1544537020393182E-01, ...
    0.4895887118324937E-02,0.5766157576611153E-03, ...
    0.3911408899973994E-02,0.1134055149899166E-01, ...
    0.2045540821981484E-01,0.1254460489172730E-01, ...
    0.5370936278809426E-02,0.5162980339438905E-02, ...
    0.5510957955189620E-02,0.1348681146014646E-01, ...
    0.5293461261550956E-02,0.2135158963906704E-01, ...
    0.3872040956962797E-02,0.4041863730898066E-01, ...
    0.3480075985572175E-02,0.3082881220603644E-01, ...
    0.3934841673430547E-02,0.4955048643733777E-02, ...
    0.5332713420324274E-02,0.3127092112918890E-01, ...
    0.3263868066569824E-01,0.3643906126193626E-02, ...
    0.1773410002914162E-02,0.2031653637703585E-01, ...
    0.9778965531187402E-02,0.4720821658621592E-02, ...
    0.2492177802034716E-02,0.2644508244904856E-01, ...
    0.1656291578773864E-01,0.9176500057842589E-02, ...
    0.1732347887472145E-01,0.1561854636309510E-01, ...
    0.2892151505244423E-02,0.3887473612274506E-01, ...
    0.3107151229520027E-02,0.2613845455837795E-01, ...
    0.3698667210865119E-01,0.2566493481700192E-02, ...
    0.1615348439627276E-01,0.3075142974941159E-01, ...
    0.2904460597417912E-02,0.2398826862979595E-01, ...
    0.7145392996633393E-02,0.8978366822810722E-02, ...
    0.2612889620260116E-01,0.1874378199539706E-01, ...
    0.1917559308928981E-01,0.4661196658021693E-02, ...
    0.1235879177909568E-01,0.4010025698621002E-01, ...
    0.3740253247098928E-01,0.3148001303878885E-02, ...
    0.3330951974159902E-01,0.1500177718625442E-01, ...
    0.6070126982071918E-02,0.3772728641321498E-01, ...
    0.1934239152098538E-01,0.2980908782028305E-01, ...
    0.1060149985804505E-01,0.1169140727671868E-01, ...
    0.2913863637730209E-01,0.2512926081791887E-02, ...
    0.4408214653268561E-01,0.2616523021560373E-01, ...
    0.2824238355893381E-02,0.3534124084616964E-01, ...
    0.3855284013804392E-01,0.3769601120920466E-01, ...
    0.3068556194531641E-01,0.1107117178829623E-01, ...
    0.1675502525295356E-01,0.4250234458185585E-02, ...
    0.2160027234544727E-01,0.4964511475611218E-01, ...
    0.3661675186805007E-01,0.1113946092997709E-01, ...
    0.4594039996241282E-02,0.2058836038539646E-01, ...
    0.8608231530163511E-02,0.1071439010120374E-01, ...
    0.1325231032613585E-01,0.2793585524817994E-01, ...
    0.1060087465526644E-01,0.4186289109506627E-01, ...
    0.2038937560848930E-01,0.4779376144301057E-02, ...
    0.4600352201496827E-01,0.4850236286164095E-01, ...
    0.1563580434067219E-01,0.2738173661329469E-01, ...
    0.3152680760913545E-01,0.4699759394230683E-01, ...
    0.1758154363368765E-01,0.1705163560097511E-01, ...
    0.4564749915571795E-01,0.3274246294590765E-01, ...
    0.9932190596257133E-02,0.3764556185875566E-01, ...
    0.5058095830353233E-01,0.4537296603002153E-01, ...
    0.4951708484325965E-01,0.3279405659962938E-02, ...
    0.9987676851635598E-02,0.4573491590126628E-01, ...
    0.9304988468918474E-02,0.1136613167624352E-01, ...
    0.5119306902579393E-01,0.4198377464395852E-01, ...
    0.1182329164018337E-01,0.1992962132997380E-01, ...
    0.2179036392120167E-01,0.2957452226482192E-01, ...
    0.3999530391333647E-01,0.1478968751484580E-01, ...
    0.3170728715126719E-01,0.1199170942901081E-01, ...
    0.4711944799135396E-01,0.2445849079621999E-01, ...
    0.1955624669135905E-01,0.4043330460830744E-01, ...
    0.2815879734180029E-02 ];

  x(1,1:n) = xs(1:n);
  x(2,1:n) = ys(1:n);
  w(1:n) = ws(1:n);

  return
end