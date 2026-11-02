%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Routine for visualization of the previously derived inverse
%$$$$$$ solution results
%$$$$$$ To be used following HeadModel_LORETA.m
%$$$$$$ 
%$$$$$$ JCM, current as of 12/29/2017
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Start file
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set path to EEGLAB parent directory
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
str = which('eeglab');
[fpath fname ext] = fileparts(str);
eeglabFolder = fpath; 

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Choose a template model for inverse solution. This should be the
%$$$$$$ same used for calculating the incerse solution
%$$$$$$
%$$$$$$ Atlas_Desikan-Killiany [488]_Colin27_2003.mat                  
%$$$$$$ Atlas_Desikan-Killiany [488]_Colin27_5003.mat                  
%$$$$$$ head_modelColin27_10003_Standard-10-5-Cap339-Destrieux148.mat  
%$$$$$$ head_modelColin27_2003_Biosemi_Standard-Cap256.mat             
%$$$$$$ head_modelColin27_2003_Standard-10-5-Cap339.mat                
%$$$$$$ head_modelColin27_2003_Standard-10-5-Cap65.mat                 
%$$$$$$ head_modelColin27_5003_Standard-10-5-Cap339-Destrieux148.mat   
%$$$$$$ head_modelColin27_5003_Standard-10-5-Cap339-ROI488.mat         
%$$$$$$ head_modelColin27_5003_Standard-10-5-Cap339.mat                
%$$$$$$ head_modelColin27_8003_Standard-10-5-Cap339.mat                
%$$$$$$ head_modelColin27_8003_Standard-10-5-Cap60.mat   
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
hmodel = 'head_modelColin27_8003_Standard-10-5-Cap339.mat';
template = [eeglabdir filesep 'plugins' filesep 'dsi' filesep 'headModel' filesep 'resources' filesep hmodel];

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Use pop_erpimage to visualize ROI source data
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Create a new EEG structure with source data in EEG.data
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
EEG_src = moveSource2DataField(EEG);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Load the headModel object as identified in the EEG structure
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
hm = headModel.loadFromFile(EEG.etc.src.hmfile);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Look into the atlas and find a region of interest
%$$$$$$
%$$$$$$ These can be found in EEG.etc.src.roi (there are 68 of them in the 
%$$$$$$ current atlas)
%$$$$$$
%$$$$$$ They are also listed at the bottom of this file
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
region = 'superiorparietal L';
ROI = find(ismember(hm.atlas.label,region));

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ First, plot single-trial source estimates
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
fig = figure;
pop_erpimage(EEG_src, 1, ROI,[],hm.atlas.label{ROI},10,1,{},[],'latency', 'yerplabel','nA/mm^2','erp','on','cbar','on');

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Next, plot the cortical map
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Compute the source ERP (collapsed within ROIs)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ERP_src_roi = mean(EEG.etc.src.act,3);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Compute the source ERP on the whole cortex
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ERP_src_cortex = mean(EEG.etc.src.actFull, 3);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Find the temporal indices surrounding event onset
%$$$$$$
%$$$$$$ Create a time index over which to compute an average
%$$$$$$
%$$$$$$ Set interval to represent this time index (with respect to time 0, 
%$$$$$$ which is the onset of the stimulus)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
interval = [150 400];
% timeOn = find(EEG_src.times == 0);
% timeOff = timeOn + interval;
time_index = find(EEG_src.times > interval(1) & EEG_src.times < interval(2));

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Get the MAX or MIN of that time window
%$$$$$$
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[~,ind_min] = min(ERP_src_roi(ROI,time_index));
[~,ind_max] = max(ERP_src_roi(ROI,time_index));

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Get the handle to the axes for the source image
%$$$$$$ 
%$$$$$$ Uncomment and use ax = get(fig,'Children'); and ax = ax(1); if
%$$$$$$ using an older version of Matlab
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ax = fig.Children(1);
% ax = get(fig,'Children');
% ax = ax(1);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Plot the cortical map on top of the existing erp-image
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
patch(...
    'vertices', hm.cortex.vertices,...
    'faces', hm.cortex.faces,...
    'FaceVertexCData', ERP_src_cortex(1:8003,time_index(ind_min)),...
    'FaceColor','interp',...
    'FaceLighting','phong',...
    'LineStyle','none',...
    'FaceAlpha',1,...
    'SpecularColorReflectance',0,...
    'SpecularExponent',25,...
    'SpecularStrength',0.25,...
    'parent',ax);
axis tight vis3d equal
set(ax,'Position',[0.1690, 0.6, 0.1673, 0.2471]);
view(ax, [-90 90])
camlight(0,180)
camlight(0,0)
title(ax,['t=' num2str(EEG_src.times(time_index(ind_min))) 'ms'])


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Interactive continuous visualization of source estimates
%$$$$$$
%$$$$$$ Visualize the source estimates on a cortical surface for the whole 
%$$$$$$ epoch
%$$$$$$ 
%$$$$$$ All required variables were defined above
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Compute channel-space ERP
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ERP = mean(EEG.data,3);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set title for figure
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$;
ftitle = [];

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set autoscale 
%$$$$$$
%$$$$$$ If false, a colorbar based on the whole data will be used
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
autoscale = 'false';

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set frames per second used to make a movie (optional)
%$$$$$$
%$$$$$$ Not checked - you might need to manually edit in the file
%$$$$$$ currentSourceViewer (see below)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
fps = 30;

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set a time range for which data are to be displayed
%$$$$$$
%$$$$$$ This can be very flexible for visualization, but not for saving 
%$$$$$$ the movie - it seems to need the entire series
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
time_line = EEG_src.times;
% time_line = [size(EEG_src.times,2)/2-100:size(EEG_src.times,2)/2+550];
% time_line = EEG_src.times(1):30:EEG_src.times(end);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Do the actual plot / movie
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
hm.plotOnModel(ERP_src_cortex, ERP, ftitle, autoscale, fps, time_line);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Note - if not using my version of the code, you will need to go
%%$$$$$ into the file "currentSourceViewer.m" and replace the embedded
%$$$$$$ function rec(obj,~,~) starting at line 301 with the following code
%$$$$$$ - remember to uncomment once replaced
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%         function rec(obj,~,~)
%             [filename,filepath] = uiputfile('*.mp4','Save movie as');
%             if isnumeric(filename);return;end
%             save_in = fullfile(filepath,filename);
%             n = size(obj.sourceMagnitud,2);
%             if obj.pointer == n
%                 obj.pointer =1;
%             end
%             obj.hAxes.XTickLabel = [];
%             obj.hAxes.YTickLabel = [];
%             obj.hAxes.ZTickLabel = [];
%             axHeight = obj.hAxes.Position(4);
%             obj.hAxes.Position(4) = 0.8;
%             axis(obj.hAxes,'on')
%             title(obj.hAxes,[obj.figureName '  ' sprintf('%f sec  (%i',obj.time(obj.pointer),obj.pointer) '/' obj.Nframes ')']);
%             frames(n) = struct('cdata',[],'colormap',[]);
%             start_loc = obj.pointer;
%             frames(obj.pointer) = getframe(obj.hFigure);
%             xx = 1;
%             while obj.pointer < n
%                 obj.next();
%                 title(obj.hAxes,[obj.figureName '  ' sprintf('%f sec  (%i',obj.time(obj.pointer),obj.pointer) '/' obj.Nframes ')']);
%                 frames(obj.pointer) = getframe(obj.hFigure);
%                 %frames2(xx,:) = frames;
%                 %xx = xx+1;
%                 pause(1/obj.fps);
%             end
%             frames(1:start_loc-1) = [];
%             %assignin('base','frames',frames2);
%             if isempty(frames), return;end
%             disp(['Now saving movie in' save_in]);
%             %movie2avi(frames, save_in, 'compression', 'None');
%             %VideoWriter(frames, save_in, 'compression', 'None');
%             v = VideoWriter(save_in,'MPEG-4');
%             v.FrameRate = 30;
%             %v.FileFormat = 'mp4';
%             open(v);
%             writeVideo(v, frames);
%             close(v);
%             title(obj.hAxes,'');
%             obj.hAxes.Position(4) = axHeight;
%             axis(obj.hAxes,'off')
%             disp('Done.')
%         end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Regions of interest available in the AAL atlas
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%     {'bankssts L'                }
%     {'bankssts R'                }
%     {'caudalanteriorcingulate L' }
%     {'caudalanteriorcingulate R' }
%     {'caudalmiddlefrontal L'     }
%     {'caudalmiddlefrontal R'     }
%     {'cuneus L'                  }
%     {'cuneus R'                  }
%     {'entorhinal L'              }
%     {'entorhinal R'              }
%     {'frontalpole L'             }
%     {'frontalpole R'             }
%     {'fusiform L'                }
%     {'fusiform R'                }
%     {'inferiorparietal L'        }
%     {'inferiorparietal R'        }
%     {'inferiortemporal L'        }
%     {'inferiortemporal R'        }
%     {'insula L'                  }
%     {'insula R'                  }
%     {'isthmuscingulate L'        }
%     {'isthmuscingulate R'        }
%     {'lateraloccipital L'        }
%     {'lateraloccipital R'        }
%     {'lateralorbitofrontal L'    }
%     {'lateralorbitofrontal R'    }
%     {'lingual L'                 }
%     {'lingual R'                 }
%     {'medialorbitofrontal L'     }
%     {'medialorbitofrontal R'     }
%     {'middletemporal L'          }
%     {'middletemporal R'          }
%     {'paracentral L'             }
%     {'paracentral R'             }
%     {'parahippocampal L'         }
%     {'parahippocampal R'         }
%     {'parsopercularis L'         }
%     {'parsopercularis R'         }
%     {'parsorbitalis L'           }
%     {'parsorbitalis R'           }
%     {'parstriangularis L'        }
%     {'parstriangularis R'        }
%     {'pericalcarine L'           }
%     {'pericalcarine R'           }
%     {'postcentral L'             }
%     {'postcentral R'             }
%     {'posteriorcingulate L'      }
%     {'posteriorcingulate R'      }
%     {'precentral L'              }
%     {'precentral R'              }
%     {'precuneus L'               }
%     {'precuneus R'               }
%     {'rostralanteriorcingulate L'}
%     {'rostralanteriorcingulate R'}
%     {'rostralmiddlefrontal L'    }
%     {'rostralmiddlefrontal R'    }
%     {'superiorfrontal L'         }
%     {'superiorfrontal R'         }
%     {'superiorparietal L'        }
%     {'superiorparietal R'        }
%     {'superiortemporal L'        }
%     {'superiortemporal R'        }
%     {'supramarginal L'           }
%     {'supramarginal R'           }
%     {'temporalpole L'            }
%     {'temporalpole R'            }
%     {'transversetemporal L'      }
%     {'transversetemporal R'      }
%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ End of file
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$