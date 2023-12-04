% testing BuildScheme and SimulateProstateSignal

%% Define VERDICT scheme

V0 = [1,2,0];
V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
V5 = [18.9, 38.8, 3000];

% Addition of b=0
Vs = [...
    V0; V1;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    V0; V5
    ];

scheme = BuildScheme(Vs);



%% simulate signal

signals = SimulateProstateSignal(scheme, fIC = 0.5, fEES = 0.5);


[bias, variance] = EvaluatefICfitting(scheme, 'RDI', fIC = 0.3, fEES = 0.5, NoiseSigma = 0.1, DEES = 2 );

