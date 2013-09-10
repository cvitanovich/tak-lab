% load knowles header file
uiopen('~/data/calib/*.mat');

% load raw data
FilterSpec='/home/andrew/data/calib/*.txt';
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec);
M=dlmread([PathName FileName],'\t'); % read tab delimited data to matrix

attens=unique(M(:,1));
% remove data with really small rms values
[rows,cols]=find(M(:,3)>=KNOWLES.coefs(3));
M=M(rows,:);
dBs = (1/KNOWLES.coefs(2))*log( (M(:,3) - KNOWLES.coefs(3)) ./ KNOWLES.coefs(1) );
M=[M dBs];
X=[];
h=figure; dim=ceil(sqrt(length(attens))); hold on;

for k=1:length(attens)
    [rows,cols]=find(M(:,1)==attens(k));
    W=M(rows,:);
    
    % select cutoffs and exclude data near noise floor or speaker limits
    [lo, hi] = select_cutoffs(log(W(:,2)),W(:,4));
    if(lo<log(32760))
        [rows,cols]=find(log(W(:,2))>lo); 
        W=W(rows,:);
        if(hi<log(32760))
            [rows,cols]=find(log(W(:,2))<hi);
            W=W(rows,:);
        else
            % not excluding any on the high end
        end
        % regression plot
        xes=log(W(:,2)); yes=W(:,4);
        figure(h); subplot(dim,dim,k); title(num2str(attens(k))); hold on;
        [RSQ, COEFFS] = regress_stats(xes,yes,0.05,[log(100) log(32760)]);
    else
        W=[]; % excluding all data at this atten
    end
    X=[X; W];
end

% attens=X(:,1);
% dbs=X(:,4);
% adj_dbs=dbs+attens;
% X=[X adj_dbs];
% figure; scatter(log(X(:,2)),adj_dbs);