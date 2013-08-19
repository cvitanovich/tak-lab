function [ILDnet] = dynsetup(num_freqs,num_ICCls_units)
%Function to set up a neural network to solve the ILDAlone specint problem
%Specifically set up for Sequential Training of a Dynamical Network (see manual)
%Architecture:
%  5-Layer network                      
%	Layers:			1		  3    5      7         9
%			Left:		Ear --> NA   VLVp   ICCls --> ICx
%						2       4  X 6  | X 8         10
%			Right:	Ear --> NA 	 VLVp   ICCls --> ICx

ILDnet                             = network; %define raw custom network
ILDnet.numLayers                   = 10; %5 layers X 2 sides (left,right)
ILDnet.numInputs						  = 2; %Left and Right ears
ILDnet.inputs{1}.size				  = num_freqs;
ILDnet.inputs{2}.size				  = num_freqs;
ILDnet.layerConnect                = [
   0 0 0 0 0 0 0 0 0 0
   0 0 0 0 0 0 0 0 0 0
   1 0 0 0 0 0 0 0 0 0
   0 1 0 0 0 0 0 0 0 0
   0 0 0 1 0 1 0 0 0 0
   0 0 1 0 1 0 0 0 0 0
   0 0 0 0 0 1 1 0 0 0
   0 0 0 0 1 0 0 1 0 0
   0 0 0 0 0 0 1 0 0 0
   0 0 0 0 0 0 0 1 0 0
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

%Layer 3 -- Left NA
layer = 3;
ILDnet.layers{layer}.size						= 20;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 4 -- Right NA
layer = 4;
ILDnet.layers{layer}.size              	= 20; 
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 5 -- Left VLVp
layer = 5;
ILDnet.layers{layer}.size						= 20;
ILDnet.layers{layer}.transferFcn       	= 'tansig';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 6 -- Right VLVp
layer = 6;
ILDnet.layers{layer}.size              	= 20;
ILDnet.layers{layer}.transferFcn       	= 'tansig';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 7 -- Left ICCls
layer = 7;
ILDnet.layers{layer}.size						= num_ICCls_units;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 8 -- Right ICCls
layer = 8;
ILDnet.layers{layer}.size              	= num_ICCls_units;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';
ILDnet.layers{layer}.topologyFcn        	= 'gridtop';

%Layer 9 -- Left ICx
layer = 9;
ILDnet.layers{layer}.size						= 1;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';

%Layer 10 -- Right ICx
layer = 10;
ILDnet.layers{layer}.size              	= 1;
ILDnet.layers{layer}.transferFcn       	= 'purelin';
ILDnet.layers{layer}.initFcn           	= 'initnw';

%Input/Output/Bias/Target
ILDnet.inputConnect(1,1)		   		 	= 1;
ILDnet.inputConnect(2,2)		     			= 1;

ILDnet.outputConnect               = [0 0 0 0 0 0 0 0 0 1];
ILDnet.biasConnect                 = [0;0;0;0;0;0;1;1;0;0];
ILDnet.targetConnect               = [0 0 0 0 0 0 0 0 0 1];

%Functions
ILDnet.initFcn                     = 'initlay';
ILDnet.performFcn                  = 'mse';
ILDnet.trainFcn                    = 'trainrp';
ILDnet.adaptFcn                    = 'adaptwb';

%Weight specifications
ILDnet.layerWeights{3,1}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{4,2}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{6,3}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{5,4}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{6,5}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{6,5}.delays			= 1; %0: do not train; 1: train
ILDnet.layerWeights{8,5}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{5,6}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{5,6}.delays			= 1; %0: do not train; 1: train
ILDnet.layerWeights{7,6}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{7,7}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{7,7}.delays    		= 1; %Recurrent connections must have nonzero delay
ILDnet.layerWeights{9,7}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{8,8}.learn			= 1; %0: do not train; 1: train
ILDnet.layerWeights{8,8}.delays    		= 1; %Recurrent connections must have nonzero delay
ILDnet.layerWeights{10,8}.learn			= 1; %0: do not train; 1: train

ILDnet.biases{7}.learn						= 1; %0: do not train; 1: train
ILDnet.biases{7}.learn						= 1; %0: do not train; 1: train
ILDnet.biases{8}.learn						= 1; %0: do not train; 1: train


%Training parameters
ILDnet.trainParam.epochs           = 500;
ILDnet.trainParam.goal             = 0.01;
ILDnet.trainParam.show             = 5;
ILDnet.trainParam.min_grad         = 1e-10;