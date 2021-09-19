% About Dataset
% images 10001 to 10348 are training images
% 10001-10202 (only Clutter) , 10203-10348 (Target + Clutter)
% images 20001 to 20348 are test images
% 20001-20146 (Target + clutter), 20147-20348 (only clutter)
%

clc
clear all
close all


  cd Database 
  
   a=0;       % variabes
   b=0;
   
   
   
  DF=[]     % create an empty array(DF Database Features)
  
  % Extraction of features from training images
  
  for i = 10001 : 10348   % number of training images
      
      i
      
      str1=int2str(i)   
      
      str2=strcat(str1,'.jpg');
      
      nor=imread(str2);    % Reading the image
      % %  colour
      
      cmap = rgb2hsv(nor); % converting RGB to HSV
      
      H=cmap(:,:,1);
      
      S=cmap(:,:,2);
      
      V=cmap(:,:,3);
      
      Hmean=mean(mean(H)); % mean Value
      Hst=std2(H);         % Standard Deviation
      
      Smean=mean(mean(S));
      Sst=std2(S);
      
      Vmean=mean(mean(V));
      Vst=std2(V);
      
      Hsk = sum(skewness(H));  % Skewness
      Ssk = sum(skewness(S));
      Vsk = sum(skewness(V));
      
      Hmin=min(imhist(H));     % Histogram minimum
      Hmax=max(imhist(H));     % Histogram maximum
      
      Smin=min(imhist(S));
      Smax=max(imhist(S));
      
      Vmin=min(imhist(V));
      Vmax=max(imhist(V));
      
      % % texture
      I=rgb2gray(nor);
      
      glcm =graycomatrix(I,'Offset',[2 0;0 2])      % Gray level coherent matrix
      
      stats1 = graycoprops(glcm,{'contrast','homogeneity'});
      
      stats2 = graycoprops(glcm,{'correlation','energy'});
      
      t1=stats1.Contrast;      % features of texture
      t2=stats1.Homogeneity;
      t3=stats2.Correlation;
      t4=stats2.Energy;
      
      % % %  shape
      % Combination of all features
      FEAT=horzcat(1,[Hmean Hst Smean Sst Vmean Vst Hsk Ssk Vsk Hmin Hmax Smin Smax Vmin Vmax t1 t2 t3 t4]);
      
      DF=[DF;FEAT];
      
      
  end
  
   
  
%   % extracting features from test image

   

 for i = 20001 : 20348  % number of test images
      
      i
      
      str1=int2str(i)   
      
      str2=strcat(str1,'.jpg');
      nor=imread(str2); 
     
   % Reading the image
      % %  colour
      
      cmap = rgb2hsv(nor); % converting RGB to HSV
      
      H=cmap(:,:,1);
      
      S=cmap(:,:,2);
      
      V=cmap(:,:,3);
      
      Hmean=mean(mean(H)); % mean Value
      Hst=std2(H);         % Standard Deviation
      
      Smean=mean(mean(S));
      Sst=std2(S);
      
      Vmean=mean(mean(V));
      Vst=std2(V);
      
      Hsk = sum(skewness(H));  % Skewness
      Ssk = sum(skewness(S));
      Vsk = sum(skewness(V));
      
      Hmin=min(imhist(H));     % Histogram minimum
      Hmax=max(imhist(H));     % Histogram maximum
      
      Smin=min(imhist(S));
      Smax=max(imhist(S));
      
      Vmin=min(imhist(V));
      Vmax=max(imhist(V));
      
      % % texture
      I=rgb2gray(nor);
      
      glcm =graycomatrix(I,'Offset',[2 0;0 2])      % Gray level coherent matrix
      
      stats1 = graycoprops(glcm,{'contrast','homogeneity'});
      
      stats2 = graycoprops(glcm,{'correlation','energy'});
      
      t1=stats1.Contrast;      % features of texture
      t2=stats1.Homogeneity;
      t3=stats2.Correlation;
      t4=stats2.Energy;
      
      % % %  shape
      % Combination of all features
      QF=horzcat(1,[Hmean Hst Smean Sst Vmean Vst Hsk Ssk Vsk Hmin Hmax Smin Smax Vmin Vmax t1 t2 t3 t4]);
      
    
      
      % % multi knn
      
      
      TrainingSet=[DF(1,:);DF(2,:);DF(3,:);DF(4,:);DF(5,:);DF(6,:);DF(7,:);DF(8,:);DF(9,:);DF(10,:); DF(11,:);DF(12,:);DF(13,:);DF(14,:);DF(15,:);DF(16,:);DF(17,:);DF(18,:);DF(19,:);DF(20,:);DF(21,:);DF(22,:);DF(23,:);DF(24,:);DF(25,:);DF(26,:);DF(27,:);DF(28,:);DF(29,:);DF(30,:);DF(31,:);DF(32,:);DF(33,:);DF(34,:);DF(35,:);DF(36,:);DF(37,:);DF(38,:);DF(39,:);DF(40,:);DF(41,:);DF(42,:);DF(43,:);DF(44,:);DF(45,:);DF(46,:);DF(47,:);DF(48,:);DF(49,:);DF(50,:);DF(51,:);DF(52,:);DF(53,:);DF(54,:);DF(55,:);DF(56,:);DF(57,:);DF(58,:);DF(59,:);DF(60,:);DF(61,:);DF(62,:);DF(63,:);DF(64,:);DF(65,:);DF(66,:);DF(67,:);DF(68,:);DF(69,:);DF(70,:);DF(71,:);DF(72,:);DF(73,:);DF(74,:);DF(75,:);DF(76,:);DF(77,:);DF(78,:);DF(79,:);DF(80,:);DF(81,:);DF(82,:);DF(83,:);DF(84,:);DF(85,:);DF(86,:);DF(87,:);DF(88,:);DF(89,:);DF(90,:);DF(91,:);DF(92,:);DF(93,:);DF(94,:);DF(95,:);DF(96,:);DF(97,:);DF(98,:);DF(99,:);DF(100,:);DF(101,:);DF(102,:);DF(103,:);DF(104,:);DF(105,:);DF(106,:);DF(107,:);DF(108,:);DF(109,:);DF(110,:); DF(111,:);DF(112,:);DF(113,:);DF(114,:);DF(115,:);DF(116,:);DF(117,:);DF(118,:);DF(119,:);DF(120,:);DF(121,:);DF(122,:);DF(123,:);DF(124,:);DF(125,:);DF(126,:);DF(127,:);DF(128,:);DF(129,:);DF(130,:);DF(131,:);DF(132,:);DF(133,:);DF(134,:);DF(135,:);DF(136,:);DF(137,:);DF(138,:);DF(139,:);DF(140,:);DF(141,:);DF(142,:);DF(143,:);DF(144,:);DF(145,:);DF(146,:);DF(147,:);DF(148,:);DF(149,:);DF(150,:);DF(151,:);DF(152,:);DF(153,:);DF(154,:);DF(155,:);DF(156,:);DF(157,:);DF(158,:);DF(159,:);DF(160,:);DF(161,:);DF(162,:);DF(163,:);DF(164,:);DF(165,:);DF(166,:);DF(167,:);DF(168,:);DF(169,:);DF(170,:);DF(171,:);DF(172,:);DF(173,:);DF(174,:);DF(175,:);DF(176,:);DF(177,:);DF(178,:);DF(179,:);DF(180,:);DF(181,:);DF(182,:);DF(183,:);DF(184,:);DF(185,:);DF(186,:);DF(187,:);DF(188,:);DF(189,:);DF(190,:);DF(191,:);DF(192,:);DF(193,:);DF(194,:);DF(195,:);DF(196,:);DF(197,:);DF(198,:);DF(199,:);DF(200,:);DF(201,:);DF(202,:);DF(203,:);DF(204,:);DF(205,:);DF(206,:);DF(207,:);DF(208,:);DF(209,:);DF(210,:); DF(211,:);DF(212,:);DF(213,:);DF(214,:);DF(215,:);DF(216,:);DF(217,:);DF(218,:);DF(219,:);DF(220,:);DF(221,:);DF(222,:);DF(223,:);DF(224,:);DF(225,:);DF(226,:);DF(227,:);DF(228,:);DF(229,:);DF(230,:);DF(231,:);DF(232,:);DF(233,:);DF(234,:);DF(235,:);DF(236,:);DF(237,:);DF(238,:);DF(239,:);DF(240,:);DF(241,:);DF(242,:);DF(243,:);DF(244,:);DF(245,:);DF(246,:);DF(247,:);DF(248,:);DF(249,:);DF(250,:);DF(251,:);DF(252,:);DF(253,:);DF(254,:);DF(255,:);DF(256,:);DF(257,:);DF(258,:);DF(259,:);DF(260,:);DF(261,:);DF(262,:);DF(263,:);DF(264,:);DF(265,:);DF(266,:);DF(267,:);DF(268,:);DF(269,:);DF(270,:);DF(271,:);DF(272,:);DF(273,:);DF(274,:);DF(275,:);DF(276,:);DF(277,:);DF(278,:);DF(279,:);DF(280,:);DF(281,:);DF(282,:);DF(283,:);DF(284,:);DF(285,:);DF(286,:);DF(287,:);DF(288,:);DF(289,:);DF(290,:);DF(291,:);DF(292,:);DF(293,:);DF(294,:);DF(295,:);DF(296,:);DF(297,:);DF(298,:);DF(299,:);DF(300,:);DF(301,:);DF(302,:);DF(303,:);DF(304,:);DF(305,:);DF(306,:);DF(307,:);DF(308,:);DF(309,:);DF(310,:); DF(311,:);DF(312,:);DF(313,:);DF(314,:);DF(315,:);DF(316,:);DF(317,:);DF(318,:);DF(319,:);DF(320,:);DF(321,:);DF(322,:);DF(323,:);DF(324,:);DF(325,:);DF(326,:);DF(327,:);DF(328,:);DF(329,:);DF(330,:);DF(331,:);DF(332,:);DF(333,:);DF(334,:);DF(335,:);DF(336,:);DF(337,:);DF(338,:);DF(339,:);DF(340,:);DF(341,:);DF(342,:);DF(343,:);DF(344,:);DF(345,:);DF(346,:);DF(347,:);DF(348,:);]
      
      % classify training set into different classes
      
      GroupTrain={'1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '1' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2' '2'} 
      
      TestSet=QF
      
      KNNModels = cell(2,1); % number of classes
      
      Y=GroupTrain
      
      classes = unique(Y);
      rng(1);  % for reproducibility
      
      % Training 
      
      for j = 1:numel(classes)
          indx = strcmp(Y',classes(j)); % crate binary classes for each classifier
           KNNModels{j} =fitcknn(DF,indx,'NumNeighbors',3,'distance','euclidean');
           
      end
      
      % Testing
      
      xGrid=QF;
      
      for j = 1:numel(classes)
          [~,score] = predict(KNNModels{j},xGrid);
          Scores(:,j) = score(:,2); % 2nd column contains positive class scores
      end


      
      [~,maxScore]= max(Scores,[],2);
      
      result=maxScore;
      
  
      
     %% counting the numbers of each class
      if i<=20146
          if result==2 % target
            a=a+1;
          end  
      end
         
       if  i>=20147 
            if  result==1 % clutter
              b=b+1;
              
              end
         
      
       end
      disp('confusion Matrix');
      disp(a);     %  correctly predicted target 
      disp(b);     %  correctly predicted clutter
      e=146-a;     %  actually target, but predictd as clutter
      f=202-b;     %  actually clutter, but predicted as target
      g=(a+b)/348;
      disp(e);
      disp(f);
      disp('Accuracy:');
      disp(g);
      disp('Probability of Detection:');
      disp((a/146));
      disp('Constant False Alarm Rate:');
      disp((f/202));
      
      cd ..
      end
      
      
     
      
      
      