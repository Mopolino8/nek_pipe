$(OBJDIR)/chkpoint.o             :src_usr/chkpoint.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/IO_tools.o             :src_usr/IO_tools.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/uparam.o               :src_usr/uparam.f;                            $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/randgen.o              :src_usr/randgen.f;                           $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/misc_stuff.o           :src_usr/misc_stuff.f;                        $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/nonlinadj.o            :src_usr/nonlinadj.f;                         $(F77) -c $(FL2) $< -o $@
