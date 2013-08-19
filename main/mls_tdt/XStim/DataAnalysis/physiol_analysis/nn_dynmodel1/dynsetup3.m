function [ILDnet] = dynsetup(num_freqs,num_ilds,num_ICCls_units)
%Function to set up a neural network to solve the ILDAlone specint problem
%Specifically set up for Batch Training of a Dynamical Network (see manual)
%Architecture:
%  2-Layer network                      
%		Input --> ICCls --> ICx --> Output

ILDnet                             = network; %define raw custom network
ILDnet.numLayers                   = 2;
ILDnet.numInputs						  = 1;
ILDnet.inputs{1}.size				  = num_freqs * num_ilds;
ILDnet.layerConnect                = [
   1 0
   1 0
];

%Layer 1 -- ICCls
layer = 1;
ILDnet.layers{layer}.size						= num_freqs * num_ilds;
ILDnet.layers{layer}.transferFcn       	= 'tansig';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.inputWeights{layer}.learn       	= 0; %0: Do NOT change, 1: CHANGE weights during training

%Layer 2 -- ICx
layer = 2;
ILDnet.layers{layer}.size              	= num_freqs;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';


%Functions
ILDnet.initFcn                     = 'initlay';
ILDnet.performFcn                  = 'mse';
ILDnet.trainFcn                    = 'trainrp';
ILDnet.adaptFcn                    = 'adaptwb';

%Weight specifications
ILDnet.layerWeights{1,1}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{2,1}.learn			= 1; %0: do not train; 1: train

ILDnet.biases{1}.learn						= 0; %0: do not train; 1: train
ILDnet.biases{2}.learn						= 0; %0: do not train; 1: train


%Training parameters
ILDnet.trainParam.epochs           = 500;
ILDnet.trainParam.goal             = 0.01;
ILDnet.trainParam.show             = 5;
ILDnet.trainParam.min_grad         = 1e-10;