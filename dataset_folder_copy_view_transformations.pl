#!/usr/bin/perl

# BigStitcher dataset copy view transformations
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# this script opens the BigStitcher dataset.xml in the working directory, and for each timepoint assuming two view setups per angle, finds the view setup with the most number of transforms for each angle (i.e. the angle used to register that view), and copies them to the other view so that view is now registered
# usage: perl dataset_folder_copy_view_transformations.pl

use Cwd qw( cwd );
use File::Spec;


my $path = Cwd::cwd();

#==================
#Main subroutine
#==================
sub main {
	my $filename = shift;
	my $html = "";
	unless ( defined($filename) && -e File::Spec->catpath( $path, $filename ) ) {
		$filename = "dataset.xml";
	}
	
	if ( open(FILE, "<" . File::Spec->catpath( $path, $filename ) ) ) {
	
		flock(FILE, LOCK_EX);
		$_ = $/;
		$/ = undef;
		$html = <FILE>; 
		flock(FILE, LOCK_UN);
		close(FILE);
		$/ = $_;
	} else {
		print "Error opening " . File::Spec->catpath( $path, $filename ) . "!\n";
		return;
	}
	my $view_setups_xml = "";
	my $copy_mode = 1; #0 means use view_setups XML content to pick angles and copy all transformations across viewsetups of the same angle but different channel/illumination settings, 1 means just copy early half of view setups to later half
	my @matched_view_setups = ();
	#print $html . "\n\n"; exit 0;
	if ( $html =~ /<ViewSetups>(.*?)<\/ViewSetups>/is ) {
		$view_setups_xml = $1;
		#$view_setups_xml =~ s/\n//g;
	}
	if ( $view_setups_xml eq "" ) {
		#$copy_mode = 1; #--already default
	} else { #$default here
		my @view_setups_list = ();
		my @view_setups = ();
		my $this_viewsetup;
		my $this_angle;
		
		while( $view_setups_xml =~ /<ViewSetup>(.*?)<\/ViewSetup>/igs ) {
			#print "While pushing\n";
			push( @view_setups, $1 );
		}
		for ( my $i=0; $i<scalar(@view_setups); $i++ ) {
			$this_viewsetup = "";
			$this_angle = "";
			if ( $view_setups[$i] =~ /<id>(.*?)<\/id>/is ) {
				$this_viewsetup = $1;
			}
			if ( $view_setups[$i] =~ /<angle>(.*?)<\/angle>/is ) {
				$this_angle = $1;
			}
			unless ( $this_viewsetup eq "" || $this_angle eq "" ) {
				#print "Pushing $this_viewsetup $this_angle\n";
				push( @view_setups_list, [$this_angle,$this_viewsetup] );
			}
		}
		
		@view_setups_list = sort{ $a->[0] cmp $b->[0] } @view_setups_list;
		my @view_setups_angles = map { $_->[0] } @view_setups_list;
		@view_setups_angles = uniq( @view_setups_angles );
		#print "View setups angles: " . join(',', @view_setups_angles) . "\n";
		for ( my $i=0; $i<scalar(@view_setups_angles); $i++ ) {
			#my @this = ( $view_setups_angles[$i] );
			my @this;
			map{ push(@this,$_->[1]) if ($_->[0] eq $view_setups_angles[$i]); } @view_setups_list;
			push( @matched_view_setups, \@this );
		}
		
		$copy_mode = 0 if ( scalar(@view_setups_list) > 0 );
	}

	print "These are the matched view setups: ";
	map{ print join(',',@{$_}) . "|"; } @matched_view_setups; 
	print "\n";
	
	#go ahead and parse XML text lines to derive the actual view setups we will ultimately copy and replace
	my @master_list; #array of arrays, 0=timepoint, 1=viewsetup, 2=line start, 3=line end
	#my $this_text;
	my @split_this_text;
	my $this_timepoint, $this_viewsetup;
	my $newline;
	my $this_start; my $this_stop;
	
	#identify and store all viewsetups for each timepoint
	while( $html =~ /<ViewRegistration\s(.*?)>(.*?)<\/ViewRegistration>/igs ) {
			$this_text = $1;
			$this_start = $-[2];
			$this_stop = $+[2];
			$this_timepoint = "";
			$this_viewsetup = "";
			
			if ( $this_stop < $this_start ) {
				print "Problem with $this_stop < $this_start, $this_text\n";
				exit 0;
			}
			
			#$this_text =~ s/\"//g;
			#print $this_text;
			if ( $this_text =~ /timepoint=\"(\d+)\"/i ) {
				$this_timepoint = $1;
			}
			if ( $this_text =~ /setup=\"(\d+)\"/i ) {
				$this_viewsetup = $1;
			}
			
			next if ( $this_viewsetup eq "" || $this_timepoint eq "" );
			#print "Pushing $this_timepoint,$this_viewsetup,$this_start,$this_stop\n";
			push( @master_list, [$this_timepoint,$this_viewsetup,$this_start,$this_stop] ) if ( defined($this_timepoint) && $this_timepoint >= 0 && defined($this_viewsetup) && $this_viewsetup >= 0 );
	}
	#exit 0;
	my @timepoint_list = map{ $_->[0] } @master_list;
	@timepoint_list = uniq(@timepoint_list);
	@timepoint_list = sort{ $a <=> $b } @timepoint_list;
	
	#check each timepoint for number of viewsetups
	my @viewsetup_list = ();
	my $add_number = 0;
	my $iteration = 0;
	for ( my $i=0; $i<scalar(@timepoint_list); $i++ ) {
		#print "Timepoint $timepoint_list[$i]...\n";
		if ( $copy_mode == 0 ) { #preferred option
			my @slave_list_lines = ();
			my $max_registration_data_line_number;
			my $max_registration_data_characters;
			my $going;
			my $newstr;
			my $add_this;
			my $added_characters;
			my $removed_characters;
			#my $newstr_characters;
			for ( my $j=0; $j<scalar(@matched_view_setups); $j++ ) {
				@slave_list_lines = ();
				$max_registration_data_characters = 0;
				$max_registration_data_line_number = -1;
				my $start_adjust_positions;
				for ( my $k=0; $k<scalar(@master_list); $k++ ) {
					if ( $master_list[$k][0] eq $timepoint_list[$i] ) {
						#right timepoint, so now look in the view setups available for this timpoint
						$going = 0;
						map{ $going = 1 if ( $_ eq $master_list[$k][1] ) } @{$matched_view_setups[$j]};
						
						next unless ( $going > 0 ); #only consider lines belonging to this matched view setups group
						#print "Checking timepoint $master_list[$k][0]\n";
						
						if ( $master_list[$k][3] - $master_list[$k][2] > $max_registration_data_characters ) { #new master viewsetup for this timepoint
							$max_registration_data_characters = $master_list[$k][3] - $master_list[$k][2];
							$max_registration_data_line_number = $k;
							push( @slave_list_lines, $k ) if ($max_registration_data_line_number>=0); #previous master relegated now to slave position
						} else {
							push( @slave_list_lines, $k );
						}
					}
				}
				#print "Looking at timepoint $timepoint_list[$i], master view setup $matched_view_setups[$j][0] / $max_registration_data_line_number, with slave lines: " . join(',',@slave_list_lines) . "\n";
				if ( $max_registration_data_line_number >= 0 && $max_registration_data_characters > 20 ) { #we're ready to copy string data
					#print "Made it here\n";
					$newstr = substr($html,$master_list[$max_registration_data_line_number][2],$master_list[$max_registration_data_line_number][3]-$master_list[$max_registration_data_line_number][2]);
					#print "Replacement text: $newstr\n\n";

					for ( my $k=0; $k<scalar(@slave_list_lines); $k++ ) {
						#first, change setup number to reflect current setup
						$add_this = $newstr;
						#$add_this =~ s/setup=\"(\d+)\"/setup=\"$slave_list_lines[$k]\"/i;
						$added_characters = length($add_this);
						
						#second, figure out where to start adjusting character offsets in $html
						$start_adjust_positions = $master_list[$slave_list_lines[$k]][3];
						
						#third, go ahead and do the replacement
						$removed_characters = $master_list[$slave_list_lines[$k]][3] - $master_list[$slave_list_lines[$k]][2];
						

						#print "Substitute after this point: " . substr($html,$start_adjust_positions,108) . "\n\n\n";
						$added_characters -= $removed_characters;
						map{ if( $_->[3] > $start_adjust_positions ) {  $_->[3] += $added_characters; }  if( $_->[2] > $start_adjust_positions ) { $_->[2] += $added_characters; } } @master_list;
						
						substr($html,$master_list[$slave_list_lines[$k]][2],$removed_characters,$add_this);

					}
					
				}
			}
		} else {
			print "Not using view setups XML content is not available at this time, aborting.\n";
			exit 0;
		
			#get viewsetups
			@viewsetup_list = ();
			map{ push( @viewsetup_list, $_->[1] ) if ( $_->[0] eq $timepoint_list[$i] ); } @master_list;
			@viewsetup_list = uniq(@viewsetup_list);
			@viewsetup_list = sort{ $a cmp $b } @viewsetup_list;
			
			if ( scalar(@viewsetup_list) % 2 == 0 ) { #even number, copy_default
				$add_number = scalar(@viewsetup_list) / 2;
				
				#copy transformations
			
			
			} else { #odd number
				print "Odd number of view setups(" . scalar(@viewsetup_list) . ") for timepoint " . $timepoint_list[$i] . ", will not copy transformations.\n";
				#$add_number = scalar(@viewsetup_list)+1) / 2;
			}
		}

	}

	rename( File::Spec->catpath( $path, $filename ), File::Spec->catpath( $path, $filename . "." . time() ) );
	if ( open(FILE, ">" . File::Spec->catpath( $path, $filename ) ) ) {
		flock(FILE, LOCK_EX);
		print FILE $html;
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error writing to " . File::Spec->catpath( $path, $filename ) . "!\n";
		return;
	}
}

#==================
#Array remove duplicates
#==================
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


#==================
#Comparison subroutines
#==================
sub max {
	my $max = shift;
	$max = $max > $_ ? $max : $_ for @_;
	return $max
}

sub min {
	my $min = shift;
	$min = $min < $_ ? $min : $_ for @_;
	return $min
}


main(@ARGV);


