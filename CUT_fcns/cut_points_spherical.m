function [xi,wi]=cut_points_spherical(N)
if N==4
    xi=[0.682259126853684,0,0;
        0,0.682259126853684,0;
        0,0,0.682259126853684;
        -0.682259126853684,0,0;
        0,-0.682259126853684,0;
        0,0,-0.682259126853684;
        0.608204882319474,0.608204882319474,0.608204882319474;
        0.608204882319474,0.608204882319474,-0.608204882319474;
        0.608204882319474,-0.608204882319474,-0.608204882319474;
        0.608204882319474,-0.608204882319474,0.608204882319474;
        -0.608204882319474,-0.608204882319474,-0.608204882319474;
        -0.608204882319474,-0.608204882319474,0.608204882319474;
        -0.608204882319474,0.608204882319474,0.608204882319474;
        -0.608204882319474,0.608204882319474,-0.608204882319474];
    wi=[0.131866518188383;
        0.131866518188383;
        0.131866518188383;
        0.131866518188383;
        0.131866518188383;
        0.131866518188383;
        0.0261001113587127;
        0.0261001113587127;
        0.0261001113587127;
        0.0261001113587127;
        0.0261001113587127;
        0.0261001113587127;
        0.0261001113587127;
        0.0261001113587127];
elseif N==6
    xi=[0,0,0;
        0.832701972673104,0,0;
        0,0.832701972673104,0;
        0,0,0.832701972673104;
        -0.832701972673104,0,0;
        0,-0.832701972673104,0;
        0,0,-0.832701972673104;
        0.747632459714073,0.747632459714073,0;
        0.747632459714073,-0.747632459714073,0;
        -0.747632459714073,0.747632459714073,0;
        -0.747632459714073,-0.747632459714073,0;
        0.747632459714073,0,0.747632459714073;
        0.747632459714073,0,-0.747632459714073;
        -0.747632459714073,0,0.747632459714073;
        -0.747632459714073,0,-0.747632459714073;
        0,0.747632459714073,0.747632459714073;
        0,0.747632459714073,-0.747632459714073;
        0,-0.747632459714073,0.747632459714073;
        0,-0.747632459714073,-0.747632459714073;
        0.429454589640814,0.429454589640814,0.429454589640814;
        0.429454589640814,0.429454589640814,-0.429454589640814;
        0.429454589640814,-0.429454589640814,0.429454589640814;
        0.429454589640814,-0.429454589640814,-0.429454589640814;
        -0.429454589640814,0.429454589640814,0.429454589640814;
        -0.429454589640814,0.429454589640814,-0.429454589640814;
        -0.429454589640814,-0.429454589640814,0.429454589640814;
        -0.429454589640814,-0.429454589640814,-0.429454589640814];
    wi=[0.0992238988178393;0.0476127646704327;0.0476127646704327;0.0476127646704327;0.0476127646704327;0.0476127646704327;0.0476127646704327;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.00908925964088841;0.0632535500154607;0.0632535500154607;0.0632535500154607;0.0632535500154607;0.0632535500154607;0.0632535500154607;0.0632535500154607;0.0632535500154607];
elseif N==8
    xi=[0,0,0;0.878243071034297,0,0;0,0.878243071034297,0;0,0,0.878243071034297;-0.878243071034297,0,0;0,-0.878243071034297,0;0,0,-0.878243071034297;0.291141229884818,0.291141229884818,0.291141229884818;0.291141229884818,0.291141229884818,-0.291141229884818;0.291141229884818,-0.291141229884818,0.291141229884818;0.291141229884818,-0.291141229884818,-0.291141229884818;-0.291141229884818,0.291141229884818,0.291141229884818;-0.291141229884818,0.291141229884818,-0.291141229884818;-0.291141229884818,-0.291141229884818,0.291141229884818;-0.291141229884818,-0.291141229884818,-0.291141229884818;0.426650784665752,0.426650784665752,0;0.426650784665752,-0.426650784665752,0;-0.426650784665752,0.426650784665752,0;-0.426650784665752,-0.426650784665752,0;0.426650784665752,0,0.426650784665752;0.426650784665752,0,-0.426650784665752;-0.426650784665752,0,0.426650784665752;-0.426650784665752,0,-0.426650784665752;0,0.426650784665752,0.426650784665752;0,0.426650784665752,-0.426650784665752;0,-0.426650784665752,0.426650784665752;0,-0.426650784665752,-0.426650784665752;0.291138122705575,0.291138122705575,0.291138122705575;0.291138122705575,0.291138122705575,-0.291138122705575;0.291138122705575,-0.291138122705575,0.291138122705575;0.291138122705575,-0.291138122705575,-0.291138122705575;-0.291138122705575,0.291138122705575,0.291138122705575;-0.291138122705575,0.291138122705575,-0.291138122705575;-0.291138122705575,-0.291138122705575,0.291138122705575;-0.291138122705575,-0.291138122705575,-0.291138122705575;0.229320049075728,0.632490645334581,0.632490645334581;0.632490645334581,0.229320049075728,0.632490645334581;0.632490645334581,0.632490645334581,0.229320049075728;-0.229320049075728,0.632490645334581,0.632490645334581;-0.632490645334581,0.229320049075728,0.632490645334581;-0.632490645334581,0.632490645334581,0.229320049075728;0.229320049075728,-0.632490645334581,0.632490645334581;0.632490645334581,-0.229320049075728,0.632490645334581;0.632490645334581,-0.632490645334581,0.229320049075728;0.229320049075728,0.632490645334581,-0.632490645334581;0.632490645334581,0.229320049075728,-0.632490645334581;0.632490645334581,0.632490645334581,-0.229320049075728;-0.229320049075728,-0.632490645334581,0.632490645334581;-0.632490645334581,-0.229320049075728,0.632490645334581;-0.632490645334581,-0.632490645334581,0.229320049075728;0.229320049075728,-0.632490645334581,-0.632490645334581;0.632490645334581,-0.229320049075728,-0.632490645334581;0.632490645334581,-0.632490645334581,-0.229320049075728;-0.229320049075728,0.632490645334581,-0.632490645334581;-0.632490645334581,0.229320049075728,-0.632490645334581;-0.632490645334581,0.632490645334581,-0.229320049075728;-0.229320049075728,-0.632490645334581,-0.632490645334581;-0.632490645334581,-0.229320049075728,-0.632490645334581;-0.632490645334581,-0.632490645334581,-0.229320049075728];
    wi=[0.00281544514784307;0.0337532853901903;0.0337532853901903;0.0337532853901903;0.0337532853901903;0.0337532853901903;0.0337532853901903;0.0140652607454182;0.0140652607454182;0.0140652607454182;0.0140652607454182;0.0140652607454182;0.0140652607454182;0.0140652607454182;0.0140652607454182;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0169138883728024;0.0140579739695725;0.0140579739695725;0.0140579739695725;0.0140579739695725;0.0140579739695725;0.0140579739695725;0.0140579739695725;0.0140579739695725;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001;0.0152796823660001];
end
end