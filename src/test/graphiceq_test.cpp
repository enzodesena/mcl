/*
 MCL
 Copyright (c) 2012-24, Enzo De Sena
 All rights reserved.
 
 Authors: Enzo De Sena, enzodesena@gmail.com
 */


#include "iirfilter.h"
#include "graphiceq.h"

namespace mcl {


bool IirFilter::PeakingFilterTest() noexcept {
  
  // Peaking filter test
  Real sample_rate(48000);
  size_t num_samples(256);
  std::vector<Real> output_data(num_samples);
  std::vector<Real> input_data(num_samples);
  input_data[0] = 1.0;
  const Real fc(500);
  const Real g(0.85);
  const Real Q(0.98);
  
  PeakingFilter peak_filter(fc, g, Q, sample_rate);
  peak_filter.Filter(&input_data[0], num_samples, &output_data[0]);
  std::vector<Real> output_data_cmp({0.994760579, -0.010091169, -0.009322856, -0.008569689, -0.007833722, -0.007116796, -0.006420542, -0.005746395, -0.005095596, -0.004469202, -0.003868099, -0.003293003, -0.002744474, -0.002222922, -0.00172862, -0.001261705, -0.000822194, -0.000409989, -2.48866E-05, 0.000333416, 0.000665311, 0.00097127, 0.001251841, 0.001507639, 0.001739336, 0.001947659, 0.00213338, 0.002297311, 0.002440295, 0.002563207, 0.002666939, 0.002752403, 0.002820522, 0.002872226, 0.002908449, 0.002930122, 0.002938172, 0.002933517, 0.002917064, 0.002889706, 0.002852317, 0.002805753, 0.002750846, 0.002688407, 0.00261922, 0.002544043, 0.002463603, 0.002378602, 0.002289709, 0.002197564, 0.002102775, 0.002005917, 0.001907536, 0.001808145, 0.001708225, 0.001608226, 0.001508566, 0.001409635, 0.001311789, 0.001215358, 0.00112064, 0.001027908, 0.000937407, 0.000849354, 0.000763942, 0.00068134, 0.000601693, 0.000525123, 0.000451733, 0.000381602, 0.000314794, 0.000251352, 0.000191303, 0.000134658, 8.14145E-05, 3.15538E-05, -1.4954E-05, -5.81511E-05, -9.80902E-05, -0.000134834, -0.000168453, -0.000199028, -0.000226644, -0.000251395, -0.000273377, -0.000292694, -0.000309452, -0.000323761, -0.000335732, -0.000345479, -0.000353118, -0.000358763, -0.000362532, -0.000364539, -0.0003649, -0.000363727, -0.000361133, -0.000357228, -0.000352119, -0.000345913, -0.00033871, -0.00033061, -0.000321711, -0.000312103, -0.000301877, -0.000291118, -0.000279907, -0.000268323, -0.000256439, -0.000244325, -0.000232048, -0.00021967, -0.000207249, -0.000194839, -0.00018249, -0.000170251, -0.000158162, -0.000146265, -0.000134594, -0.000123182, -0.000112059, -0.000101249, -9.07766E-05, -8.06603E-05, -7.09173E-05, -6.15619E-05, -5.26057E-05, -4.40577E-05, -3.59248E-05, -2.82116E-05, -2.09207E-05, -1.40526E-05, -7.60624E-06, -1.5788E-06, 4.0341E-06, 9.23821E-06, 1.40406E-05, 1.84495E-05, 2.24741E-05, 2.61247E-05, 2.94123E-05, 3.23487E-05, 3.49462E-05, 3.72179E-05, 3.91771E-05, 4.08375E-05, 4.22132E-05, 4.33183E-05, 4.41672E-05, 4.47743E-05, 4.51539E-05, 4.53204E-05, 4.52879E-05, 4.50706E-05, 4.46822E-05, 4.41363E-05, 4.34462E-05, 4.26247E-05, 4.16844E-05, 4.06376E-05, 3.94959E-05, 3.82708E-05, 3.69731E-05, 3.56133E-05, 3.42013E-05, 3.27466E-05, 3.12582E-05, 2.97446E-05, 2.82138E-05, 2.66733E-05, 2.51303E-05, 2.35912E-05, 2.20621E-05, 2.05486E-05, 1.9056E-05, 1.75889E-05, 1.61516E-05, 1.47479E-05, 1.33814E-05, 1.2055E-05, 1.07715E-05, 9.53313E-06, 8.34187E-06, 7.19935E-06, 6.1069E-06, 5.06553E-06, 4.07597E-06, 3.1387E-06, 2.25394E-06, 1.42167E-06, 6.41661E-07, -8.65048E-08, -7.63445E-07, -1.38994E-06, -1.96692E-06, -2.49547E-06, -2.97678E-06, -3.41217E-06, -3.80304E-06, -4.15089E-06, -4.45728E-06, -4.72385E-06, -4.95227E-06, -5.14427E-06, -5.3016E-06, -5.42603E-06, -5.51935E-06, -5.58333E-06, -5.61977E-06, -5.63045E-06, -5.6171E-06, -5.58148E-06, -5.52528E-06, -5.45017E-06, -5.35779E-06, -5.24972E-06, -5.1275E-06, -4.99263E-06, -4.84655E-06, -4.69065E-06, -4.52625E-06, -4.35463E-06, -4.17701E-06, -3.99453E-06, -3.8083E-06, -3.61933E-06, -3.42861E-06, -3.23704E-06, -3.04548E-06, -2.85472E-06, -2.66548E-06, -2.47845E-06, -2.29424E-06, -2.11342E-06, -1.9365E-06, -1.76394E-06, -1.59614E-06, -1.43346E-06, -1.27623E-06, -1.12471E-06, -9.7912E-07, -8.39658E-07, -7.06469E-07, -5.79665E-07, -4.59324E-07, -3.45491E-07, -2.38183E-07, -1.37387E-07, -4.30651E-08, 4.48458E-08, 1.2643E-07, 2.01792E-07, 2.71056E-07, 3.34361E-07, 3.91862E-07, 4.43726E-07});
  ASSERT(IsEqual(output_data_cmp, output_data));
  
  return true;
}



bool IirFilter::PeakHighShelfTest() noexcept {
  
  // Peaking filter test
  Real sample_rate(48000);
  size_t num_samples(256);
  std::vector<Real> output_data(num_samples);
  std::vector<Real> input_data(num_samples);
  input_data[0] = 1.0;
  const Real fc(500);
  const Real g(0.85);
  const Real Q(0.98);
  
  PeakHighShelf peak_filter(fc, g, Q, sample_rate);
  peak_filter.Filter(&input_data[0], num_samples, &output_data[0]);
  std::vector<Real> output_data_cmp({0.852377516, 0.004885228, 0.005127942, 0.005335967, 0.005510661, 0.005653424, 0.005765694, 0.005848935, 0.005904631, 0.005934279, 0.005939383, 0.005921448, 0.005881973, 0.005822447, 0.005744342, 0.005649111, 0.005538182, 0.005412953, 0.005274793, 0.005125032, 0.004964964, 0.004795841, 0.004618872, 0.004435223, 0.004246009, 0.004052301, 0.003855119, 0.003655432, 0.003454161, 0.003252173, 0.003050286, 0.002849266, 0.002649829, 0.002452638, 0.002258309, 0.002067408, 0.001880451, 0.00169791, 0.001520208, 0.001347725, 0.001180797, 0.001019716, 0.000864738, 0.000716074, 0.000573903, 0.000438364, 0.000309565, 0.000187579, 7.24508E-05, -3.58057E-05, -0.000137203, -0.000231779, -0.000319595, -0.000400736, -0.000475305, -0.000543425, -0.000605236, -0.000660894, -0.000710567, -0.000754437, -0.000792697, -0.000825548, -0.000853201, -0.000875873, -0.000893786, -0.000907167, -0.000916247, -0.000921258, -0.000922433, -0.000920007, -0.000914211, -0.000905278, -0.000893436, -0.000878911, -0.000861926, -0.000842699, -0.000821441, -0.000798361, -0.000773659, -0.000747532, -0.000720167, -0.000691747, -0.000662444, -0.000632427, -0.000601855, -0.000570878, -0.000539641, -0.00050828, -0.000476921, -0.000445684, -0.000414683, -0.00038402, -0.000353792, -0.000324087, -0.000294988, -0.000266566, -0.00023889, -0.000212018, -0.000186004, -0.000160894, -0.000136727, -0.000113538, -9.13548E-05, -7.01994E-05, -5.00891E-05, -3.1036E-05, -1.30472E-05, 3.87441E-06, 1.97305E-05, 3.45264E-05, 4.82715E-05, 6.09784E-05, 7.2663E-05, 8.3344E-05, 9.30429E-05, 0.000101784, 0.000109592, 0.000116497, 0.000122526, 0.000127713, 0.000132089, 0.000135688, 0.000138545, 0.000140694, 0.000142171, 0.000143013, 0.000143255, 0.000142934, 0.000142087, 0.000140748, 0.000138954, 0.00013674, 0.00013414, 0.000131189, 0.000127919, 0.000124362, 0.000120551, 0.000116515, 0.000112284, 0.000107886, 0.000103348, 9.86972E-05, 9.39573E-05, 8.91523E-05, 8.43046E-05, 7.94354E-05, 7.45646E-05, 6.97111E-05, 6.48922E-05, 6.01243E-05, 5.54225E-05, 5.08006E-05, 4.62713E-05, 4.18463E-05, 3.7536E-05, 3.33497E-05, 2.92958E-05, 2.53815E-05, 2.16133E-05, 1.79963E-05, 1.4535E-05, 1.12331E-05, 8.09329E-06, 5.11746E-06, 2.30683E-06, -3.38086E-07, -2.81747E-06, -5.13212E-06, -7.2834E-06, -9.27324E-06, -1.1104E-05, -1.27787E-05, -1.43004E-05, -1.5673E-05, -1.69003E-05, -1.79868E-05, -1.8937E-05, -1.97558E-05, -2.04481E-05, -2.10192E-05, -2.14745E-05, -2.18193E-05, -2.20592E-05, -2.21999E-05, -2.22468E-05, -2.22058E-05, -2.20823E-05, -2.18821E-05, -2.16105E-05, -2.12732E-05, -2.08754E-05, -2.04224E-05, -1.99194E-05, -1.93715E-05, -1.87835E-05, -1.81601E-05, -1.7506E-05, -1.68255E-05, -1.61229E-05, -1.54023E-05, -1.46674E-05, -1.39221E-05, -1.31698E-05, -1.24139E-05, -1.16574E-05, -1.09032E-05, -1.01542E-05, -9.41289E-06, -8.68156E-06, -7.96244E-06, -7.2575E-06, -6.56858E-06, -5.89731E-06, -5.24516E-06, -4.61344E-06, -4.00329E-06, -3.41572E-06, -2.85157E-06, -2.31154E-06, -1.7962E-06, -1.30599E-06, -8.41223E-07, -4.02099E-07, 1.12958E-08, 3.98977E-07, 7.61058E-07, 1.09774E-06, 1.40932E-06, 1.69616E-06, 1.9587E-06, 2.19745E-06, 2.41296E-06, 2.60585E-06, 2.7768E-06, 2.92652E-06, 3.05573E-06, 3.16524E-06, 3.25583E-06, 3.32835E-06, 3.38363E-06, 3.42253E-06, 3.44592E-06, 3.45467E-06, 3.44967E-06, 3.43178E-06, 3.40187E-06, 3.36081E-06, 3.30943E-06, 3.24859E-06, 3.17909E-06, 3.10175E-06, 3.01735E-06, 2.92664E-06, 2.83037E-06, 2.72925E-06, 2.62397E-06});
  ASSERT(IsEqual(output_data_cmp, output_data));

  return true;
}


bool IirFilter::PeakLowShelfTest() noexcept {
  
  // Peaking filter test
  Real sample_rate(48000);
  size_t num_samples(256);
  std::vector<Real> output_data(num_samples);
  std::vector<Real> input_data(num_samples);
  input_data[0] = 1.0;
  const Real fc(500);
  const Real g(0.85);
  const Real Q(0.98);
  
  PeakLowShelf peak_filter(fc, g, Q, sample_rate);
  peak_filter.Filter(&input_data[0], num_samples, &output_data[0]);
  std::vector<Real> output_data_cmp({0.997210725, -0.00571531, -0.005966509, -0.006174058, -0.006339956, -0.006466255, -0.00655505, -0.00660846, -0.006628622, -0.006617679, -0.006577769, -0.006511017, -0.006419525, -0.006305366, -0.006170574, -0.00601714, -0.005847005, -0.005662055, -0.005464114, -0.005254943, -0.005036235, -0.004809613, -0.004576626, -0.004338745, -0.004097369, -0.003853815, -0.003609322, -0.003365051, -0.003122083, -0.002881421, -0.002643991, -0.002410641, -0.002182144, -0.0019592, -0.001742437, -0.001532412, -0.001329616, -0.001134473, -0.000947344, -0.000768533, -0.000598282, -0.00043678, -0.000284166, -0.000140526, -5.90305E-06, 0.000119704, 0.000236338, 0.000344077, 0.000443036, 0.000533359, 0.000615223, 0.000688826, 0.000754394, 0.000812173, 0.000862425, 0.000905432, 0.000941486, 0.000970893, 0.000993967, 0.001011031, 0.00102241, 0.001028436, 0.001029442, 0.00102576, 0.001017721, 0.001005654, 0.000989884, 0.000970731, 0.000948507, 0.000923518, 0.000896063, 0.00086643, 0.000834899, 0.000801738, 0.000767206, 0.00073155, 0.000695006, 0.000657796, 0.000620133, 0.000582217, 0.000544234, 0.00050636, 0.000468757, 0.000431575, 0.000394954, 0.00035902, 0.000323888, 0.000289662, 0.000256434, 0.000224288, 0.000193293, 0.000163513, 0.000135, 0.000107795, 8.19341E-05, 5.74422E-05, 3.43375E-05, 1.26305E-05, -7.67521E-06, -2.65827E-05, -4.41009E-05, -6.02446E-05, -7.50336E-05, -8.84925E-05, -0.00010065, -0.00011154, -0.000121197, -0.000129662, -0.000136976, -0.000143185, -0.000148333, -0.00015247, -0.000155645, -0.000157908, -0.000159311, -0.000159904, -0.000159739, -0.000158869, -0.000157344, -0.000155215, -0.000152533, -0.000149346, -0.000145703, -0.000141651, -0.000137235, -0.0001325, -0.000127488, -0.00012224, -0.000116797, -0.000111194, -0.000105469, -9.96553E-05, -9.37846E-05, -8.78872E-05, -8.19915E-05, -7.61238E-05, -7.03085E-05, -6.45683E-05, -5.89238E-05, -5.33939E-05, -4.79958E-05, -4.27449E-05, -3.76547E-05, -3.27375E-05, -2.80036E-05, -2.34619E-05, -1.912E-05, -1.49839E-05, -1.10584E-05, -7.34685E-06, -3.85164E-06, -5.73892E-07, 2.48628E-06, 5.32976E-06, 7.95833E-06, 1.03746E-05, 1.25821E-05, 1.45848E-05, 1.63875E-05, 1.79957E-05, 1.9415E-05, 2.06518E-05, 2.17129E-05, 2.26053E-05, 2.33363E-05, 2.39134E-05, 2.43445E-05, 2.46373E-05, 2.47999E-05, 2.48403E-05, 2.47664E-05, 2.45864E-05, 2.43081E-05, 2.39394E-05, 2.3488E-05, 2.29615E-05, 2.23673E-05, 2.17127E-05, 2.10045E-05, 2.02496E-05, 1.94546E-05, 1.86256E-05, 1.77688E-05, 1.68897E-05, 1.59938E-05, 1.50864E-05, 1.41721E-05, 1.32557E-05, 1.23414E-05, 1.1433E-05, 1.05344E-05, 9.64878E-06, 8.77938E-06, 7.92896E-06, 7.10008E-06, 6.29501E-06, 5.51575E-06, 4.76407E-06, 4.0415E-06, 3.34933E-06, 2.6886E-06, 2.06019E-06, 1.46474E-06, 9.02712E-07, 3.74378E-07, -1.2015E-07, -5.80925E-07, -1.00814E-06, -1.40214E-06, -1.76339E-06, -2.09245E-06, -2.39003E-06, -2.65689E-06, -2.8939E-06, -3.102E-06, -3.28221E-06, -3.43557E-06, -3.56322E-06, -3.6663E-06, -3.74601E-06, -3.80354E-06, -3.84015E-06, -3.85705E-06, -3.85551E-06, -3.83677E-06, -3.80207E-06, -3.75263E-06, -3.68967E-06, -3.61438E-06, -3.52793E-06, -3.43144E-06, -3.32604E-06, -3.21279E-06, -3.09272E-06, -2.96684E-06, -2.83609E-06, -2.7014E-06, -2.56363E-06, -2.42361E-06, -2.28211E-06, -2.13988E-06, -1.99759E-06, -1.8559E-06, -1.71539E-06, -1.57662E-06, -1.4401E-06, -1.30628E-06, -1.17558E-06, -1.04839E-06, -9.25034E-07, -8.05812E-07, -6.90981E-07, -5.80761E-07, -4.75338E-07, -3.74862E-07, -2.79452E-07});
  ASSERT(IsEqual(output_data_cmp, output_data));

  return true;
}

//bool IirFilter::GraphicEqTest() noexcept {
//  
//  // Graphic EQ test
//  
//  std::vector<Real> fc({ 250.0, 500.0, 1000.0, 2000.0, 4000.0 });
//  std::vector<Real> gains({0.9, 0.85, 0.88, 0.8, 0.75});
//  Real Q = 0.98;
//  int fs = 48e3;
//  int num_samples = 256;
//  
//  std::vector<Real> output_data(num_samples);
//  std::vector<Real> input_data(num_samples);
//  input_data[0] = 1.0;
//  
//  GraphicEq eq(gains, fc, Q, fs);
//  eq.Filter(&input_data[0], num_samples, &output_data[0]);
//
//  std::vector<Real> output_data_cmp({0.789794889,  -0.003782045, 0.010303789, 0.015598704, 0.016070502, 0.014404195, 0.011990325, 0.009409928, 0.006905777, 0.004641463, 0.002764293, 0.001371886, 0.000472902, -1.49692E-05, -0.000230256, -0.000318545, -0.000392073, -0.000510555, -0.000683466, -0.00088612, -0.001079845, -0.001228757, -0.001309534, -0.001314137, -0.001247528, -0.001123023, -0.00095745, -0.000767407, -0.000567046, -0.000367201, -0.000175465, 3.25643E-06, 0.000166057, 0.000311396, 0.000438616, 0.000547623, 0.000638684, 0.000712316, 0.000769229, 0.000810289, 0.000836493, 0.000848941, 0.000848812, 0.00083733, 0.000815732, 0.000785243, 0.000747048, 0.00070228, 0.000652002, 0.00059721, 0.000538826, 0.000477706, 0.00041464, 0.000350358, 0.000285533, 0.000220783, 0.000156673, 9.3715E-05, 3.23695E-05, -2.69575E-05, -8.39128E-05, -0.000138197, -0.000189563, -0.000237817, -0.000282813, -0.000324455, -0.000362692, -0.000397514, -0.000428949, -0.00045706, -0.00048194, -0.000503705, -0.000522496, -0.000538468, -0.000551788, -0.000562634, -0.000571189, -0.000577634, -0.000582155, -0.000584928, -0.000586128, -0.00058592, -0.000584461, -0.000581897, -0.000578364, -0.000573987, -0.000568878, -0.000563138, -0.000556859, -0.000550121, -0.000542991, -0.000535532, -0.000527792, -0.000519816, -0.000511638, -0.000503287, -0.000494786, -0.000486153, -0.000477402, -0.000468543, -0.000459583, -0.000450528, -0.000441379, -0.00043214, -0.000422812, -0.000413395, -0.000403891, -0.000394302, -0.000384628, -0.000374873, -0.00036504, -0.000355133, -0.000345158, -0.000335122, -0.000325031, -0.000314895, -0.000304722, -0.000294524, -0.000284312, -0.000274096, -0.00026389, -0.000253706, -0.000243558, -0.000233458, -0.00022342, -0.000213457, -0.000203583, -0.000193809, -0.00018415, -0.000174615, -0.000165218, -0.000155968, -0.000146876, -0.000137952, -0.000129204, -0.00012064, -0.000112267, -0.000104092, -9.61207E-05, -8.8357E-05, -8.08055E-05, -7.34694E-05, -6.63512E-05, -5.9453E-05, -5.27758E-05, -4.63203E-05, -4.00865E-05, -3.40739E-05, -2.82816E-05, -2.27079E-05, -1.73511E-05, -1.22089E-05, -7.27869E-06, -2.55762E-06, 1.95746E-06, 6.2699E-06, 1.03832E-05, 1.43012E-05, 1.80275E-05, 2.15662E-05, 2.49213E-05, 2.80968E-05, 3.10968E-05, 3.39255E-05, 3.6587E-05, 3.90853E-05, 4.14246E-05, 4.3609E-05, 4.56425E-05, 4.75291E-05, 4.92727E-05, 5.08773E-05, 5.23467E-05, 5.36847E-05, 5.48951E-05, 5.59815E-05, 5.69476E-05, 5.7797E-05, 5.85332E-05, 5.91596E-05, 5.96799E-05, 6.00972E-05, 6.0415E-05, 6.06367E-05, 6.07653E-05, 6.08042E-05, 6.07566E-05, 6.06256E-05, 6.04143E-05, 6.01258E-05, 5.97631E-05, 5.93292E-05, 5.88271E-05, 5.82598E-05, 5.76301E-05, 5.69409E-05, 5.61951E-05, 5.53954E-05, 5.45446E-05, 5.36454E-05, 5.27005E-05, 5.17126E-05, 5.06842E-05, 4.96179E-05, 4.85163E-05, 4.73818E-05, 4.62169E-05, 4.50239E-05, 4.38053E-05, 4.25633E-05, 4.13001E-05, 4.00181E-05, 3.87193E-05, 3.74058E-05, 3.60798E-05, 3.47433E-05, 3.33981E-05, 3.20463E-05, 3.06897E-05, 2.93302E-05, 2.79693E-05, 2.6609E-05, 2.52508E-05, 2.38964E-05, 2.25473E-05, 2.12049E-05, 1.98709E-05, 1.85465E-05, 1.72331E-05, 1.59321E-05, 1.46446E-05, 1.33719E-05, 1.2115E-05, 1.08752E-05, 9.65344E-06, 8.45071E-06, 7.26798E-06, 6.10613E-06, 4.96603E-06, 3.84847E-06, 2.75421E-06, 1.68396E-06, 6.38376E-07, -3.81922E-07, -1.37636E-06, -2.34443E-06, -3.28563E-06, -4.19953E-06, -5.08573E-06, -5.94389E-06, -6.77369E-06, -7.57487E-06, -8.34718E-06, -9.09045E-06, -9.80452E-06, -1.04893E-05});
//  ASSERT(IsEqual(output_data_cmp, output_data));
//  ASSERT(false);
//  return true;
//}

bool IirFilter::GraphicEqTest() noexcept {
  
  // Graphic EQ test
  
  std::vector<Real> fc({ 250.0, 500.0, 1000.0, 2000.0, 4000.0 });
  std::vector<Real> gains({0.9, 0.85, 0.88, 0.8, 0.75});
  Real Q = 0.98;
  Real fs = 48e3;
  size_t num_samples = 256;
  
  std::vector<Real> output_data(num_samples);
  std::vector<Real> input_data(num_samples);
  input_data[0] = 1.0;
  
  GraphicEq eq(gains, fc, Q, fs);
  eq.Filter(&input_data[0], num_samples, &output_data[0]);

  std::vector<Real> output_data_cmp({0.789794889,  -0.003782045, 0.010303789, 0.015598704, 0.016070502, 0.014404195, 0.011990325, 0.009409928, 0.006905777, 0.004641463, 0.002764293, 0.001371886, 0.000472902, -1.49692E-05, -0.000230256, -0.000318545, -0.000392073, -0.000510555, -0.000683466, -0.00088612, -0.001079845, -0.001228757, -0.001309534, -0.001314137, -0.001247528, -0.001123023, -0.00095745, -0.000767407, -0.000567046, -0.000367201, -0.000175465, 3.25643E-06, 0.000166057, 0.000311396, 0.000438616, 0.000547623, 0.000638684, 0.000712316, 0.000769229, 0.000810289, 0.000836493, 0.000848941, 0.000848812, 0.00083733, 0.000815732, 0.000785243, 0.000747048, 0.00070228, 0.000652002, 0.00059721, 0.000538826, 0.000477706, 0.00041464, 0.000350358, 0.000285533, 0.000220783, 0.000156673, 9.3715E-05, 3.23695E-05, -2.69575E-05, -8.39128E-05, -0.000138197, -0.000189563, -0.000237817, -0.000282813, -0.000324455, -0.000362692, -0.000397514, -0.000428949, -0.00045706, -0.00048194, -0.000503705, -0.000522496, -0.000538468, -0.000551788, -0.000562634, -0.000571189, -0.000577634, -0.000582155, -0.000584928, -0.000586128, -0.00058592, -0.000584461, -0.000581897, -0.000578364, -0.000573987, -0.000568878, -0.000563138, -0.000556859, -0.000550121, -0.000542991, -0.000535532, -0.000527792, -0.000519816, -0.000511638, -0.000503287, -0.000494786, -0.000486153, -0.000477402, -0.000468543, -0.000459583, -0.000450528, -0.000441379, -0.00043214, -0.000422812, -0.000413395, -0.000403891, -0.000394302, -0.000384628, -0.000374873, -0.00036504, -0.000355133, -0.000345158, -0.000335122, -0.000325031, -0.000314895, -0.000304722, -0.000294524, -0.000284312, -0.000274096, -0.00026389, -0.000253706, -0.000243558, -0.000233458, -0.00022342, -0.000213457, -0.000203583, -0.000193809, -0.00018415, -0.000174615, -0.000165218, -0.000155968, -0.000146876, -0.000137952, -0.000129204, -0.00012064, -0.000112267, -0.000104092, -9.61207E-05, -8.8357E-05, -8.08055E-05, -7.34694E-05, -6.63512E-05, -5.9453E-05, -5.27758E-05, -4.63203E-05, -4.00865E-05, -3.40739E-05, -2.82816E-05, -2.27079E-05, -1.73511E-05, -1.22089E-05, -7.27869E-06, -2.55762E-06, 1.95746E-06, 6.2699E-06, 1.03832E-05, 1.43012E-05, 1.80275E-05, 2.15662E-05, 2.49213E-05, 2.80968E-05, 3.10968E-05, 3.39255E-05, 3.6587E-05, 3.90853E-05, 4.14246E-05, 4.3609E-05, 4.56425E-05, 4.75291E-05, 4.92727E-05, 5.08773E-05, 5.23467E-05, 5.36847E-05, 5.48951E-05, 5.59815E-05, 5.69476E-05, 5.7797E-05, 5.85332E-05, 5.91596E-05, 5.96799E-05, 6.00972E-05, 6.0415E-05, 6.06367E-05, 6.07653E-05, 6.08042E-05, 6.07566E-05, 6.06256E-05, 6.04143E-05, 6.01258E-05, 5.97631E-05, 5.93292E-05, 5.88271E-05, 5.82598E-05, 5.76301E-05, 5.69409E-05, 5.61951E-05, 5.53954E-05, 5.45446E-05, 5.36454E-05, 5.27005E-05, 5.17126E-05, 5.06842E-05, 4.96179E-05, 4.85163E-05, 4.73818E-05, 4.62169E-05, 4.50239E-05, 4.38053E-05, 4.25633E-05, 4.13001E-05, 4.00181E-05, 3.87193E-05, 3.74058E-05, 3.60798E-05, 3.47433E-05, 3.33981E-05, 3.20463E-05, 3.06897E-05, 2.93302E-05, 2.79693E-05, 2.6609E-05, 2.52508E-05, 2.38964E-05, 2.25473E-05, 2.12049E-05, 1.98709E-05, 1.85465E-05, 1.72331E-05, 1.59321E-05, 1.46446E-05, 1.33719E-05, 1.2115E-05, 1.08752E-05, 9.65344E-06, 8.45071E-06, 7.26798E-06, 6.10613E-06, 4.96603E-06, 3.84847E-06, 2.75421E-06, 1.68396E-06, 6.38376E-07, -3.81922E-07, -1.37636E-06, -2.34443E-06, -3.28563E-06, -4.19953E-06, -5.08573E-06, -5.94389E-06, -6.77369E-06, -7.57487E-06, -8.34718E-06, -9.09045E-06, -9.80452E-06, -1.04893E-05});
  ASSERT(IsEqual(output_data_cmp, output_data)); // TODO: for some reason the serialised filter does not decay as fast as the "normal" one. Double check.
  
  return true;
}

} // namespace mcl
