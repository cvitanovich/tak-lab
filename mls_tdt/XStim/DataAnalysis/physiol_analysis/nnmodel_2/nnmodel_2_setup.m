%Script to set up a neural network to solve the ILDAlone specint problem
%Specifically set up for Batch Training of a Static Network (see manual)
%Architecture:
%  2-Layer network
%    Input(16x21 ILD/freq spec) --> Layer 1 (numfreqs X numILDs) --> Layer 2 (1 ICCls or ICx unit)
%                               \ /                              \ /
%                         Weights are 1 if               Weights are a function of
%                         ILD/freq matches,              the Euclidean distance of a location's
%                         0 otherwise                    ILD spectrum from the unit's optimal
%                                                        ILD spectrum

clear;
ICspacenet                             = network; %define raw custom network
ICspacenet.numLayers                   = 2;
ICspacenet.layerConnect                = [0 0
                                          1 0];

%Layer 1
ICspacenet.layers{1}.transferFcn       = 'tansig';
ICspacenet.layers{1}.initFcn           = '';
ICspacenet.layers{1}.topologyFcn       = 'gridtop';

%Layer 2
ICspacenet.layers{2}.size              = 1; %Layer representing the ICCls/ICx unit
ICspacenet.layers{2}.transferFcn       = 'purelin';
ICspacenet.layers{2}.initFcn           = '';

%Input/Output/Bias/Target
ICspacenet.numInputs                   = 1;
ICspacenet.inputConnect(1,:)           = 1;
ICspacenet.inputConnect(2,:)           = 0;
ICspacenet.outputConnect               = [0 1];
ICspacenet.biasConnect                 = [0;0];
ICspacenet.targetConnect               = [0 1];

%Functions
ICspacenet.initFcn                     = 'initlay';
ICspacenet.performFcn                  = 'mse';
ICspacenet.trainFcn                    = 'trainrp';
ICspacenet.adaptFcn                    = 'adaptwb';

%Weight specifications
ICspacenet.inputWeights{1}.learn       = 1; %0: Do NOT change, 1: CHANGE weights during training
ICspacenet.inputWeights{1}.learnFcn    = '';
ICspacenet.layerWeights{1,1}.learn     = 1; %0: Do NOT change, 1: CHANGE weights during training
ICspacenet.layerWeights{1,1}.delays    = 1; %Recurrent connections must have nonzero delay
ICspacenet.layerWeights{1,1}.learnFcn  = '';
ICspacenet.layerWeights{2,1}.learn     = 0; %0: Do NOT change, 1: CHANGE weights during training
ICspacenet.layerWeights{2,1}.learnFcn  = '';

%Training parameters
ICspacenet.trainParam.epochs           = 500;
ICspacenet.trainParam.goal             = 0.01;
ICspacenet.trainParam.show             = 5;
ICspacenet.trainParam.min_grad         = 1e-10;

%Save the network
save d:\mlspezio\matlab\save\ICspacenet_init ICspacenet
