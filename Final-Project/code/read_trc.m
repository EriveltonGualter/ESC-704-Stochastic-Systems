function [ trc ] = read_trc(raw)
% trc = read_trc(raw)   
% 
%     This function creates the structure data for MotionAnalysis data
% 
%     INPUTS:
%         raw = [n, 86] = Raw Data 
%
%     OUTPUTS:
%         trc = Tracker motion
%             .frm = [n, 1] = frame
%             .t = [n, 1] = tempo [s]
%             .IJ = Incisura jugular
%             .PX = Processo xifoide 
%             .C7 = Sétima vertebra cervical
%             .T8 = Oitava vertebra torácica 
%             .RAC/.LAC = Articulação acrômio clavicular direita/esquerda
%             .RAS/.RPS/.RAI/.RPI = Cluster braço direito
%             .LAS/.LPS/.LAI/.LPI = Cluster braço esquerdo
%             .RLE/.LLE = Epicôndilo lateral do úmero direita/esquerda
%             .RME/.LME = Epicôndilo medial do úmero direita/esquerda
%             .RRS/.LRS = Processo estiloide do rádio  
%             .RUS/.LUS = Processo estiloide da ulna 
%             .R3M/.L3M = Cabeça do terceiro osso metacarpo direita/esquerda
%             .R5M/.L5M = Cabeça do quinto osso metacarpo  
%             .RW/.LW = Eixo roda traseira direita/esquerda

    for i=1:length(raw)
        trc.frm(i,1) = raw(i,1);
        trc.t(i,1) = raw(i,2);

        trc.IJ.x(i,1) = raw(i,3);
        trc.IJ.y(i,2) = raw(i,4);
        trc.IJ.z(i,3) = raw(i,5);

        trc.PX.x(i,1) = raw(i,6);
        trc.PX.y(i,2) = raw(i,7);
        trc.PX.z(i,3) = raw(i,8);

        trc.C7.x(i,1) = raw(i,9);
        trc.C7.y(i,2) = raw(i,10);
        trc.C7.z(i,3) = raw(i,11);

        trc.T8.x(i,1) = raw(i,12);
        trc.T8.y(i,2) = raw(i,13);
        trc.T8.z(i,3) = raw(i,14);

        trc.RAC.x(i,1) = raw(i,15);
        trc.RAC.y(i,2) = raw(i,16);
        trc.RAC.z(i,3) = raw(i,17);

        trc.RAS.x(i,1) = raw(i,18);
        trc.RAS.y(i,2) = raw(i,19);
        trc.RAS.z(i,3) = raw(i,20);

        trc.RPS.x(i,1) = raw(i,21);
        trc.RPS.y(i,2) = raw(i,22);
        trc.RPS.z(i,3) = raw(i,23);

        trc.RAI.x(i,1) = raw(i,24);
        trc.RAI.y(i,2) = raw(i,25);
        trc.RAI.z(i,3) = raw(i,26);

        trc.RPI.x(i,1) = raw(i,27);
        trc.RPI.y(i,2) = raw(i,28);
        trc.RPI.z(i,3) = raw(i,29);

        trc.RLE.x(i,1) = raw(i,30);
        trc.RLE.y(i,2) = raw(i,31);
        trc.RLE.z(i,3) = raw(i,32);

        trc.RME.x(i,1) = raw(i,33);
        trc.RME.y(i,2) = raw(i,34);
        trc.RME.z(i,3) = raw(i,35);

        trc.RRS.x(i,1) = raw(i,36);
        trc.RRS.y(i,2) = raw(i,37);
        trc.RRS.z(i,3) = raw(i,38);

        trc.RUS.x(i,1) = raw(i,39);
        trc.RUS.y(i,2) = raw(i,40);
        trc.RUS.z(i,3) = raw(i,41);

        trc.R3M.x(i,1) = raw(i,42);
        trc.R3M.y(i,2) = raw(i,43);
        trc.R3M.z(i,3) = raw(i,44);

        trc.R5M.x(i,1) = raw(i,45);
        trc.R5M.y(i,2) = raw(i,46);
        trc.R5M.z(i,3) = raw(i,47);

        trc.LAC.x(i,1) = raw(i,48);
        trc.LAC.y(i,2) = raw(i,49);
        trc.LAC.z(i,3) = raw(i,50);

        trc.LAS.x(i,1) = raw(i,51);
        trc.LAS.y(i,2) = raw(i,52);
        trc.LAS.z(i,3) = raw(i,53);

        trc.LPS.x(i,1) = raw(i,54);
        trc.LPS.y(i,2) = raw(i,55);
        trc.LPS.z(i,3) = raw(i,56);

        trc.LAI.x(i,1) = raw(i,57);
        trc.LAI.y(i,2) = raw(i,58);
        trc.LAI.z(i,3) = raw(i,59);

        trc.LPI.x(i,1) = raw(i,60);
        trc.LPI.y(i,2) = raw(i,61);
        trc.LPI.z(i,3) = raw(i,62);

        trc.LLE.x(i,1) = raw(i,63);
        trc.LLE.y(i,2) = raw(i,64);
        trc.LLE.z(i,3) = raw(i,65);

        trc.LME.x(i,1) = raw(i,66);
        trc.LME.y(i,2) = raw(i,67);
        trc.LME.z(i,3) = raw(i,68);

        trc.LRS.x(i,1) = raw(i,69);
        trc.LRS.y(i,2) = raw(i,70);
        trc.LRS.z(i,3) = raw(i,71);

        trc.LUS.x(i,1) = raw(i,72);
        trc.LUS.y(i,2) = raw(i,73);
        trc.LUS.z(i,3) = raw(i,74);

        trc.L3M.x(i,1) = raw(i,75);
        trc.L3M.y(i,2) = raw(i,76);
        trc.L3M.z(i,3) = raw(i,77);

        trc.L5M.x(i,1) = raw(i,78);
        trc.L5M.y(i,2) = raw(i,79);
        trc.L5M.z(i,3) = raw(i,80);

        trc.RW.x(i,1) = raw(i,81);
        trc.RW.y(i,2) = raw(i,82);
        trc.RW.z(i,3) = raw(i,83);

        trc.LW.x(i,1) = raw(i,84);
        trc.LW.y(i,2) = raw(i,85);
        trc.LW.z(i,3) = raw(i,86);
    end
end