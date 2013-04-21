function [ILDnet] = dynsetup(num_freqs,num_ICCls_units)
%Function to set up a neural network to solve the ILDAlone specint problem
%Specifically set up for Sequential Training of a Dynamical Network (see manual)
%Architecture:
%  5-Layer network                      
%	Layers:			1	      3         5
%			Left:		Ear -->  ICCls --> ICx
%						2        4			 6
%			Right:	Ear -->  ICCls --> ICx

ILDnet                             = network; %define raw custom network
ILDnet.numLayers                   = 6; %3 layers X 2 sides (left,right)
ILDnet.numInputs						  = 2; %Left and Right ears
ILDnet.inputs{1}.size				  = num_freqs;
ILDnet.inputs{2}.size				  = num_freqs;
ILDnet.layerConnect                = [
   0 0 0 0 0 0
   0 0 0 0 0 0
   1 0 1 0 0 0
   0 1 0 1 0 0
   0 0 1 0 0 0
   0 0 0 1 0 0
];

%Layer 1 -- Left Ear
layer = 1;
ILDnet.layers{layer}.size						= num_freqs;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.inputWeights{layer}.learn       	= 0; %0: Do NOT change, 1: CHANGE weights during training

%Layer 2 -- Right Ear
layer = 2;
ILDnet.layers{layer}.size              	= num_freqs;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.inputWeights{layer}.learn       	= 0; %0: Do NOT change, 1: CHANGE weights during training

%Layer 3 -- Left ICCls
layer = 3;
ILDnet.layers{layer}.size						= num_ICCls_units;
ILDnet.layers{layer}.transferFcn       	= 'tansig';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 4 -- Right ICCls
layer = 4;
ILDnet.layers{layer}.size              	= num_ICCls_units;
ILDnet.layers{layer}.transferFcn       	= 'tansig';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 5 -- Left ICx
layer = 5;
ILDnet.layers{layer}.size						= 1;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';

%Layer 6 -- Right ICx
layer = 6;
ILDnet.layers{layer}.size              	= 1;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';

%Input/Output/Bias/Target
ILDnet.inputConnect(1,1)		   		 	= 1;
ILDnet.inputConnect(2,2)		     			= 1;

ILDnet.outputConnect               = [0 0 0 0 0 1];
ILDnet.biasConnect                 = [0;0;1;1;0;0];
ILDnet.targetConnect               = [0 0 0 0 0 1];

%Functions
ILDnet.initFcn                     = 'initlay';
ILDnet.performFcn                  = 'mse';
ILDnet.trainFcn                    = 'trainrp';
ILDnet.adaptFcn                    = 'adaptwb';

%Weight specifications
ILDnet.layerWeights{3,1}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{4,2}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{3,3}.learn			= 1;
ILDnet.layerWeights{3,3}.delays			= 1;
ILDnet.layerWeights{5,3}.learn			= 1;
ILDnet.layerWeights{4,4}.learn			= 1;
ILDnet.layerWeights{4,4}.delays			= 1;
ILDnet.layerWeights{6,4}.learn			= 1;

ILDnet.biases{3}.learn						= 1; %0: do not train; 1: train
ILDnet.biases{4}.learn						= 1; %0: do not train; 1: train


%Training parameters
ILDnet.trainParam.epochs           = 500;
ILDnet.trainParam.goal             = 0.01;
ILDnet.trainParam.show             = 5;
ILDnet.trainParam.min_grad         = 1e-10;