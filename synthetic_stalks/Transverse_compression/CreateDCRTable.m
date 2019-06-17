function [Stalk_TableDCR] = CreateDCRTable(Table,range,ext_xDCR,ext_yDCR,ext_tDCR,ext_rhoDCR,int_xDCR,int_yDCR,int_tDCR,int_rhoDCR)
% CreateDCRTable.m: Take the variables created from PrepSections_V2.m and
% make a new data table, copied from SMALL_CURVES_V2_3_1500.mat

Stalk_TableDCR = Table(range(1):range(2),:);

Stalk_TableDCR = removevars(Stalk_TableDCR,{'Ext_X','Ext_Y','Int_X','Int_Y','xbar','ybar'});

N = size(Stalk_TableDCR,1);

ext_xDCR = single(ext_xDCR);
ext_yDCR = single(ext_yDCR);
ext_tDCR = single(ext_tDCR);
ext_rhoDCR = single(ext_rhoDCR);
int_xDCR = single(int_xDCR);
int_yDCR = single(int_yDCR);
int_tDCR = single(int_tDCR);
int_rhoDCR = single(int_rhoDCR);

Ext_X = cell(N,1);
Ext_Y = cell(N,1);
Ext_T = cell(N,1);
Ext_Rho = cell(N,1);
Int_X = cell(N,1);
Int_Y = cell(N,1);
Int_T = cell(N,1);
Int_Rho = cell(N,1);

for i = 1:N
    Ext_X{i} = ext_xDCR(:,i);
    Ext_Y{i} = ext_yDCR(:,i);
    Ext_T{i} = ext_tDCR(:,i);
    Ext_Rho{i} = ext_rhoDCR(:,i);
    Int_X{i} = int_xDCR(:,i);
    Int_Y{i} = int_yDCR(:,i);
    Int_T{i} = int_tDCR(:,i);
    Int_Rho{i} = int_rhoDCR(:,i);   
end


% Save new boundary profiles into Stalk_TableDCR
Stalk_TableDCR = addvars(Stalk_TableDCR,Ext_X,'Before','rind_t');
Stalk_TableDCR = addvars(Stalk_TableDCR,Ext_Y,'After','Ext_X');
Stalk_TableDCR = addvars(Stalk_TableDCR,Int_X,'After','Ext_Y');
Stalk_TableDCR = addvars(Stalk_TableDCR,Int_Y,'After','Int_X');
Stalk_TableDCR = addvars(Stalk_TableDCR,Ext_T,'After','Int_Y');
Stalk_TableDCR = addvars(Stalk_TableDCR,Ext_Rho,'After','Ext_T');
Stalk_TableDCR = addvars(Stalk_TableDCR,Int_T,'After','Ext_Rho');
Stalk_TableDCR = addvars(Stalk_TableDCR,Int_Rho,'After','Int_T');

end