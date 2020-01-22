clc;clear all; close all;
Qs = 1:7;
QsStr = cellstr([repmat('qusetion ',length(Qs),1) num2str(Qs')]);
while 1
    pdfOrsol = menu("Open solution or questions :","Questions (PDF)","Solutions");
    clc;close all;
    if pdfOrsol == 0
        disp("Ok! Good bye :)");
        break
    elseif pdfOrsol == 1
        disp("Questions file is openning...")
        system(fullfile("..","BSP_1_98_2m.pdf"));
        continue
    end
    while 1
        chosen = menu("Choose the question :",QsStr);
        clc;close all;
        if chosen == 0
            disp("Ok! Good bye :)");
            break
        else
            run(['Q' num2str(Qs(chosen)) '/main.m'])
        end
    end
end