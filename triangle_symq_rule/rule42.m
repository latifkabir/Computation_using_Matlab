function [ x, y, w ] = rule42 ( )

%*****************************************************************************80
%
%% RULE42 returns the rule of degree 42.
%
%  Discussion:
%
%    Order 42 (324 pts)
%    1/6 data for 42-th order quadrature with 62 nodes.
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
       0.00000000000000000000000000000000, ...
       0.00000000000000000000000000000000, ...
      -0.33732186702983072527377866348931, ...
       0.00000000000000000000000000000000, ...
      -0.96033929602597311557761977477877, ...
      -0.86369982899641926289375115271667, ...
       0.00000000000000000000000000000000, ...
       0.00000000000000000000000000000000, ...
      -0.91964505433817813657264837251525E-01, ...
      -0.66309355251219477483888932315586, ...
      -0.51372916997472296488281957020097, ...
       0.00000000000000000000000000000000, ...
      -0.93438263780905817601029667101954, ...
       0.00000000000000000000000000000000, ...
      -0.90377238041680626119957295786590, ...
      -0.97636560454756210194148963029962, ...
      -0.80976763015790175270655294347859, ...
      -0.11503943729553541704016281376984, ...
      -0.45642869772220040454862543963416, ...
      -0.26081632436550845472752474357332, ...
      -0.23737442933070842287301049364327, ...
      -0.48002697917900107156245025485608, ...
      -0.55650162501787387537981069378810, ...
      -0.12478176574417714212774268424213, ...
       0.00000000000000000000000000000000, ...
      -0.57250293774060918335092055220851, ...
      -0.35082904042609762405060413462916, ...
       0.00000000000000000000000000000000, ...
       0.00000000000000000000000000000000, ...
      -0.48329385573259627109084340701262, ...
      -0.60723425879311620189016063272093, ...
      -0.74479369220077260375364623234338, ...
      -0.88899121520171718881192636679144E-01, ...
      -0.64093927012580804843536900337703, ...
      -0.72840614316212545054184733409861, ...
      -0.34632334976132410062989430769368, ...
      -0.29328359092418193020774881585849, ...
      -0.52075640841362497935400311743224, ...
      -0.60555027935120159340828726711903, ...
      -0.40905454399991287433613498814895, ...
       0.00000000000000000000000000000000, ...
      -0.68384877689179067855591079553230, ...
      -0.71845211882128619528987877814782, ...
      -0.21667354415688649871182254929595, ...
      -0.78401417474182939603497549800273, ...
      -0.73829866676407222514323372659105E-01, ...
      -0.13454529284893910979363670735151, ...
      -0.81224627731204864480024020247250, ...
      -0.21779874120415795023777409358666, ...
       0.00000000000000000000000000000000, ...
      -0.42173742894104082835679644928203, ...
      -0.86724517518683599270071051238882, ...
      -0.35974948772649811066870719069685, ...
      -0.38462016485447150927117414922040, ...
      -0.14709579616219243877526898806227, ...
      -0.14236698354849056304093448822305, ...
       0.00000000000000000000000000000000, ...
      -0.21806859696894763488778637237838, ...
      -0.28058519641950796678974392628735, ...
       0.00000000000000000000000000000000, ...
       0.00000000000000000000000000000000, ...
       0.00000000000000000000000000000000 ];
  y = [ ... ...
      -0.23915669727967103118380139768186, ...
       0.11499706797515916998220806401861E+01, ...
      -0.57463021225838519900335995927075, ...
      -0.56879992979770302521884041359830, ...
      -0.54196210567272946906489712244628, ...
      -0.53941219280079277870643403803930, ...
      -0.50772285124871792168226530442572, ...
       0.92445442054483368113851672487764, ...
      -0.57476336247113177371961348404936, ...
      -0.48132277490581309659947493263723, ...
      -0.37176697642006607510778844628504, ...
       0.10049252920684353614918795242733E+01, ...
      -0.57493742331160895521040050048537, ...
      -0.54376585322681251262386416118604, ...
      -0.55951105488338182442194601355457, ...
      -0.57287342016015178837509330570583, ...
      -0.50897875397137854857739836711645, ...
      -0.55605671624505328368102486115945, ...
      -0.57286458799044539871933065364668, ...
      -0.51230069551399328623142182353309, ...
      -0.54774164573700082169258700885382, ...
      -0.42881250819630884111523888550059, ...
      -0.47691618600212933305894604508553, ...
      -0.51872351401291148614495359408056, ...
       0.10611738733054864308352869822745E+01, ...
      -0.57297391133539792503817569661015, ...
      -0.55889259681253444754878375605606, ...
       0.67490581305685926714955055515407, ...
      -0.46308066751251361359457711328632, ...
      -0.55408939280498580695400859722506, ...
      -0.42433782075470535918840983649065, ...
      -0.52032748775510677492310958477448, ...
      -0.23994545417022376699165559718980, ...
      -0.52285819108589385499701905311092, ...
      -0.46066493316460819615558728243469, ...
      -0.25420978730877827119837106988602, ...
      -0.46670083015487351222116455572352, ...
      -0.52093525817617790107842375316544, ...
      -0.55435992266257116866227670316197, ...
      -0.34138644830088341585231904393856, ...
       0.78549511301950868330086537396862, ...
      -0.57282030241995550160692247503933, ...
      -0.55400012643472741059711642382680, ...
      -0.24499725128675631298290375529787, ...
      -0.57293824798152066163974973440469, ...
      -0.40234854190615793941968454639826, ...
      -0.14476929626336181753272157091150, ...
      -0.55359930548468954074315532763614, ...
      -0.57155725996820959476202397056243, ...
       0.74817466669102061182393603775266E-01, ...
      -0.48091898059119731232834971500976, ...
      -0.57249401685345627316687746770457, ...
      -0.41046630772117739922581443734300, ...
      -0.52677593302205091340352395053193, ...
      -0.46723591737438723999494490231575, ...
      -0.32805865585485880306160348707314, ...
       0.30647862202829743938620359376295, ...
      -0.40551576066406470095805909037152, ...
      -0.33418952446785743728949888770460, ...
       0.54905382739186092457188995702651, ...
      -0.32644308263379735102228951317930, ...
      -0.14180931586841959084250934615779 ];
  w = [ ... ...
       0.23336887717611856261356883070280E-02, ...
       0.22846739164921267903695204331260E-04, ...
       0.70111433879772498743791597802441E-03, ...
       0.66736115582873727562875088212462E-03, ...
       0.57785786514264065438532593664687E-03, ...
       0.13648017145261238754376858784678E-02, ...
       0.18715796130141111886024261596391E-02, ...
       0.12311750904528985305433912829374E-02, ...
       0.73641178059201489942713406653212E-03, ...
       0.33400912254975143418927736967276E-02, ...
       0.53408360185582866807409529152528E-02, ...
       0.87691656602387380470263190564516E-03, ...
       0.29827840345921013265121144680296E-03, ...
       0.14428656899060768897862657143696E-02, ...
       0.89736946555330205057357879483966E-03, ...
       0.24622832560853176787327086769666E-03, ...
       0.22197294779296558696560399072253E-02, ...
       0.25038963454662955321871600907872E-02, ...
       0.10178082444714138863786539439139E-02, ...
       0.41750972622189718039985200715331E-02, ...
       0.28539390170943817441698264694008E-02, ...
       0.53618401327742042860547051137744E-02, ...
       0.44379187212093563379576658435829E-02, ...
       0.43923515022689111826190793996648E-02, ...
       0.62671836400915488143041776852213E-03, ...
       0.97955449417192964190566775890366E-03, ...
       0.23018208943548581359392383454856E-02, ...
       0.28689426702522553567921914933795E-02, ...
       0.29825739568587811806411249714095E-02, ...
       0.25047104675084795836734599265446E-02, ...
       0.52092111037958286914269299684586E-02, ...
       0.28471024101248386546695067189237E-02, ...
       0.82474954978823005019031886199947E-02, ...
       0.31999957236382169770094914684218E-02, ...
       0.37199859798352839350442391912941E-02, ...
       0.85220480529686853126175155374059E-02, ...
       0.55913173485803385482484428528564E-02, ...
       0.38298884403760618569504568998379E-02, ...
       0.23089127874100000149196844935704E-02, ...
       0.70117980296802485289704436077025E-02, ...
       0.24321512738535918412817989755737E-02, ...
       0.93669661545711090880311728905074E-03, ...
       0.20642879082447440946641719056348E-02, ...
       0.90936551597467743554493293539138E-02, ...
       0.79986620899126243521988101233587E-03, ...
       0.77120700714742987371071118791067E-02, ...
       0.10219559553814186449346549825429E-01, ...
       0.17241106961996584893369591793848E-02, ...
       0.13912182187991904328859295409100E-02, ...
       0.53283478006725665670848896066726E-02, ...
       0.53185152093575276335966863060532E-02, ...
       0.68994902614553096607967297139181E-03, ...
       0.68499915239657437963246941954227E-02, ...
       0.39920622096982370641881855622790E-02, ...
       0.64073752422181504921770295563876E-02, ...
       0.86898690635024727222144915241917E-02, ...
       0.47974852921883757762301535128063E-02, ...
       0.74399462871562674711612756662098E-02, ...
       0.81490368631924462501235690710386E-02, ...
       0.39847607162807112945480365181103E-02, ...
       0.44476356973374278433219749730315E-02, ...
       0.52129984983515280965046257082396E-02 ];

  return
end
