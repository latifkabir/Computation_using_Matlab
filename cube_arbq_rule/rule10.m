function [ x, w ] = rule10 ( n )

%*****************************************************************************80
%
%% RULE10 returns the rule of degree 10.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 July 2014
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
%    Output, real X(3,N), the coordinates of the nodes.
%
%    Output, real W(N), the weights.
%
  xs = [ ...
    -0.9290008391594218E+00, 0.4559324934613976E+00, ...
    -0.8503323372224584E+00, -0.6258718000763578E+00, ...
    -0.7657521094877171E+00, 0.5194124679448002E+00, ...
    -0.1876523482951094E-01, 0.7781249769033519E+00, ...
    0.7907391104275550E+00, 0.1931186493871086E-01, ...
    0.7973052188120806E-01, 0.8940978313791650E+00, ...
    -0.3969913823357744E+00, -0.5840731119485252E+00, ...
    -0.3054403489443687E+00, -0.6250814570086429E+00, ...
    -0.4273170491243774E+00, 0.8887351719726477E+00, ...
    -0.1425902095867877E+00, 0.1007385849466660E+00, ...
    0.9914273856893784E+00, -0.2110597230830005E+00, ...
    -0.9344242484094948E+00, 0.9254796648599809E+00, ...
    0.5871962886124195E+00, -0.5644486965116591E+00, ...
    -0.3450006491762153E+00, -0.1925706173126046E+00, ...
    0.6150095319350374E+00, -0.9410296227058319E+00, ...
    0.3990765859738513E+00, 0.3675511374298872E+00, ...
    0.9468259750893930E+00, -0.8323347102578033E+00, ...
    0.9487564158646586E+00, 0.6918585765278729E-01, ...
    0.7257633299128992E+00, -0.3310667096656268E+00, ...
    0.3501561292254360E+00, -0.2984756141933171E+00, ...
    0.7448197310529086E+00, 0.8904795041635960E-01, ...
    -0.8524054154380974E+00, -0.9797468974234023E+00, ...
    0.2069917638147098E+00, -0.9222673238089781E+00, ...
    -0.7404783544230075E+00, 0.4462428682423911E+00, ...
    0.2499725669820603E+00, -0.5978592395842010E+00, ...
    0.3880147299726153E+00, 0.9540597896735835E+00, ...
    -0.1204215373750881E+00, 0.5635716786856889E+00, ...
    -0.6198158466921444E+00, 0.5423781905790268E+00, ...
    0.2041052298129538E+00, 0.8615414935964518E+00, ...
    0.9625611083095826E+00, -0.5668586728612423E+00, ...
    -0.9565835991972781E+00, 0.9290333570079982E+00, ...
    0.9603761783535766E+00, -0.9001174228174761E+00, ...
    0.7409829822994354E+00, -0.9186328686279674E+00, ...
    0.7988074968670636E+00, -0.8055374206467150E+00, ...
    -0.2303785439930381E+00, -0.7165954822802608E+00, ...
    0.7098003466268242E+00, -0.2370105920082338E+00, ...
    -0.9973208321714092E+00 ];
  ys = [ ...
    0.5674482687326375E+00, -0.2911223108663519E+00, ...
    0.1975635105828826E+00, -0.6882273878784136E+00, ...
    -0.9313977287710445E+00, -0.5339817326366849E+00, ...
    0.5593088502110057E+00, 0.6938671668589624E+00, ...
    0.7457324218656869E+00, -0.4591331058869646E-01, ...
    0.5440689882793791E+00, -0.8616026786491190E+00, ...
    0.7739289299329076E+00, -0.2114642239345504E+00, ...
    -0.3131293885912573E-01, 0.9186402247572539E+00, ...
    0.9877228633852757E+00, -0.8848784715526166E+00, ...
    -0.8430339428373445E+00, -0.9121169131818918E+00, ...
    0.8352108151249428E-01, -0.5108468353410851E+00, ...
    0.7934672273226390E+00, 0.5843397847611630E-02, ...
    0.1982224584099894E+00, 0.6529491925552273E+00, ...
    0.3453783893258016E+00, -0.6965527925846071E+00, ...
    -0.6030971553224019E+00, -0.9773971815452341E+00, ...
    0.9306828192128160E+00, 0.1729189273247773E+00, ...
    0.9423175295395478E+00, 0.3117471500716507E+00, ...
    0.8484099142032112E+00, -0.4549656193618034E+00, ...
    0.4364950089820075E+00, -0.2119737649763384E-01, ...
    -0.7501696781614142E-01, 0.3377085041383336E+00, ...
    0.5119386721405460E+00, -0.2861289818382340E+00, ...
    -0.9435758438429453E+00, -0.5501358135716665E+00, ...
    -0.8151761180652929E+00, -0.4382794463039796E+00, ...
    0.1037779030310166E+00, -0.9327537847317822E+00, ...
    0.7781924335471115E+00, -0.6696667136365475E+00, ...
    0.9993628918113781E+00, 0.6532557117383904E+00, ...
    0.8321683233238897E+00, 0.8340145881278882E+00, ...
    -0.3990216364491641E+00, -0.7026952932032947E+00, ...
    0.4592299742407868E+00, -0.4282015920212867E+00, ...
    0.2004721513667550E+00, -0.8934779235340296E+00, ...
    0.8942502209608590E+00, -0.8738897670233189E+00, ...
    -0.5705526001256189E+00, -0.1379917314855049E+00, ...
    -0.1665061725494578E+00, -0.7946830985787409E+00, ...
    0.9767806321337382E+00, 0.9597550703525793E+00, ...
    0.9603148869900205E+00, 0.6723297616898366E+00, ...
    -0.9804400598708566E+00, -0.9949665182949334E+00, ...
    0.4191106171156644E+00 ];
  zs = [ ...
    -0.9933978864296309E+00, 0.9901550709677349E+00, ...
    0.9670485250396404E+00, 0.9686386295527810E+00, ...
    -0.9482984392355279E+00, -0.9781218402197481E+00, ...
    -0.9668923135113735E+00, 0.9569807092335304E+00, ...
    -0.9255203116327040E+00, -0.8549048193170219E+00, ...
    0.9206982447618198E+00, 0.9211244686937629E+00, ...
    -0.7317954338646777E+00, -0.9289922028324820E+00, ...
    0.8753848332253249E+00, -0.8882119979065144E+00, ...
    0.9878899913703325E+00, -0.8880356823987890E+00, ...
    -0.8903407082095307E+00, 0.8876467178896359E+00, ...
    -0.9463172167846347E+00, 0.7267804723276519E+00, ...
    0.8338672039848504E+00, 0.8514657617058597E+00, ...
    -0.8105375103476049E+00, 0.7553648186259019E+00, ...
    -0.5850428048791508E+00, -0.1054105129348611E+00, ...
    0.7438442367847414E+00, 0.8658386707444876E+00, ...
    0.7640316811327564E+00, 0.6612763241373982E+00, ...
    0.7191463486728814E+00, -0.7362541429366616E+00, ...
    -0.6578813706843174E+00, -0.6173927832617907E+00, ...
    0.5688448669683742E+00, -0.2424692431476031E+00, ...
    -0.3220397059827600E+00, 0.2741662573733258E+00, ...
    -0.4245710530722387E+00, 0.3779936020945108E+00, ...
    -0.4890743188433689E+00, -0.8427007062812649E+00, ...
    0.3110361053541443E+00, 0.7265683900043480E+00, ...
    0.5114587028650807E+00, -0.6590556333897338E+00, ...
    -0.5656407036165850E+00, -0.5801968526465481E+00, ...
    -0.8814901582749934E+00, 0.3049921985976026E+00, ...
    0.4604030935118299E+00, 0.1644483100006628E+00, ...
    0.1392355684867057E+00, -0.1804588400155518E+00, ...
    -0.1659257531038997E-01, -0.6408450861337968E+00, ...
    -0.2283594679036603E+00, 0.5423020975283729E+00, ...
    -0.5488666370681280E+00, -0.1935317969193018E+00, ...
    0.4390831472700056E+00, -0.3197962461756297E+00, ...
    0.1476513989767390E+00, 0.1104669057496399E+00, ...
    -0.1823430458926981E+00, 0.3668414192631035E+00, ...
    -0.2174815141367195E+00, -0.1332596398840131E+00, ...
    0.4442289032147039E+00, -0.1887850946760386E+00, ...
    0.2545601348113754E+00 ];
  ws = [ ...
    0.7213361346970690E-02, 0.2333141303565080E-01, ...
    0.1679006811176461E-01, 0.1914958033825363E-01, ...
    0.7581033822063699E-02, 0.2434386571403711E-01, ...
    0.2909344157431100E-01, 0.1700037166331166E-01, ...
    0.1673233204518952E-01, 0.3947746719767090E-01, ...
    0.3908321158531827E-01, 0.1015197562376277E-01, ...
    0.2887452492906294E-01, 0.4051167157357858E-01, ...
    0.5023577704859834E-01, 0.1585721420277247E-01, ...
    0.7185787768043105E-02, 0.1081642505043930E-01, ...
    0.3035331237698869E-01, 0.2526839495335092E-01, ...
    0.8807297819671642E-02, 0.5593986671876056E-01, ...
    0.1373773962571757E-01, 0.2661615086380097E-01, ...
    0.5999969968705832E-01, 0.4972097477068091E-01, ...
    0.6999983849508175E-01, 0.5598146788784706E-01, ...
    0.5112821739971896E-01, 0.4853200937844245E-02, ...
    0.2643720405354587E-01, 0.7488081622291798E-01, ...
    0.8050832752047557E-02, 0.4291659968333356E-01, ...
    0.1294625305397554E-01, 0.7403978691481330E-01, ...
    0.5117832177948896E-01, 0.7869167027903835E-01, ...
    0.7987922983101066E-01, 0.7935975973606037E-01, ...
    0.5190634933748844E-01, 0.8887073925421256E-01, ...
    0.1520654567592403E-01, 0.1437138357007375E-01, ...
    0.5539120615704855E-01, 0.3308981979630965E-01, ...
    0.6553371660845619E-01, 0.3169478052156868E-01, ...
    0.6276849569730045E-01, 0.6211610103098184E-01, ...
    0.1053947144115475E-01, 0.2331767126681452E-01, ...
    0.6027945648014443E-01, 0.4978984285341077E-01, ...
    0.7974451931076965E-01, 0.6688248399941439E-01, ...
    0.9777913635381968E-01, 0.5098104982036934E-01, ...
    0.2951246816860325E-01, 0.4353390845153333E-01, ...
    0.1466253257039164E-01, 0.2081467883348036E-01, ...
    0.2633254465708832E-01, 0.5423792444480278E-01, ...
    0.8745492482441419E-01, 0.3119590040598747E-01, ...
    0.1563786213616101E-01, 0.2105985290616827E-01, ...
    0.3523362981701141E-01, 0.7539771307788777E-01, ...
    0.1972567979376486E-01, 0.2513903851381411E-01, ...
    0.2400953849626596E-01 ];

  x(1,1:n) = xs(1:n);
  x(2,1:n) = ys(1:n);
  x(3,1:n) = zs(1:n);
  w(1:n) = ws(1:n);

  return
end
