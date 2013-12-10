PREFIX                  = @prefix@
BINDIR                  = @bindir@
LIBDIR                  = @libdir@
INCLUDEDIR              = @includedir@

CC                      = @cc@
CFLAGS                  = @cflags@
BASE_CFLAGS             = @base_cflags@
LDFLAGS                 = @ldflags@
OPTFLAGS                = @optflags@

AR                      = @ar@
CP                      = @cp@
RM                      = @rm@
RANLIB                  = @ranlib@
MKDIR                   = @mkdir@

#------------------------------------------------------------------------------#
COBJFILEPOSTFIX         = @cobjfilepostfix@

TARGET                  = @fstrexec_targetfile@

#------------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .c

.c.@cobjfilepostfix@:
	$(CC) -c $(BASE_CFLAGS) $(CFLAGS) $(OPTFLAGS) $< -o $@

#------------------------------------------------------------------------------#
OBJS = \
	starter.@cobjfilepostfix@

HEADERS =

#------------------------------------------------------------------------------#
all: @build_target@

build-default: $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(LDFLAGS)

$(OBJS): $(HEADERS)

install:
	if [ ! -d $(PREFIX) ]; then \
		$(MKDIR) $(PREFIX); \
	fi
	if [ ! -d $(BINDIR) ]; then \
		$(MKDIR) $(BINDIR); \
	fi
#	$(CP) ./$(TARGET) $(BINDIR)/fistr
#	$(CP) ./$(TARGET) $(BINDIR)/hecmw_part
#	$(CP) ./$(TARGET) $(BINDIR)/hecmw_vis

clean:
	$(RM) *.$(COBJFILEPOSTFIX)

distclean: clean
	$(RM) $(TARGET)
