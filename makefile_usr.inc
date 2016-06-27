$(OBJDIR)/chkpoint.o             :src_usr/chkpoint.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/IO_tools.o             :src_usr/IO_tools.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/uparam.o               :src_usr/uparam.f;                            $(F77) -c $(FL2) $< -o $@
