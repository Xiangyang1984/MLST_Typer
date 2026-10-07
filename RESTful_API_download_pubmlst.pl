#!/usr/bin/env perl
# Generalized from Pneumo-Typer/ST_tool/RESTful_API_download_pubmlst.pl.
# Keeps the original PubMLST OAuth 1.0A workflow, but accepts any PubMLST
# sequence-definition database URL and stores ONE global token set in auth/.

use strict;
use warnings;
use 5.010;
use FindBin;
use lib "$FindBin::Bin/lib";
use Net::OAuth 0.20;
$Net::OAuth::PROTOCOL_VERSION = Net::OAuth::PROTOCOL_VERSION_1_0A;
use HTTP::Request::Common;
use LWP::UserAgent;
use JSON qw(decode_json);
use Data::Random qw(rand_chars);
use Config::Tiny;
use Getopt::Long qw(:config no_ignore_case);
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Cwd qw(abs_path);

use constant CONSUMER_KEY    => 'A0IPaxPdVT84GyXeeXG8Iup2';
use constant CONSUMER_SECRET => 'fgdmRPXCAHFOaMrghi2wZe57G3JkuIfiTOOcVxP2xL';

my %opts = (
    auth_dir => "$FindBin::Bin/auth",
);
GetOptions(
    'a|arguments=s' => \$opts{'a'},
    'm|method=s'    => \$opts{'m'},
    'r|route=s'     => \$opts{'r'},
    'u|url=s'       => \$opts{'u'},
    'b|base-url=s'  => \$opts{'b'},
    'o|output=s'    => \$opts{'o'},
    'auth-dir=s'    => \$opts{'auth_dir'},
    'h|help'        => \$opts{'h'},
) or die("Error in command line arguments\n");

if ($opts{'h'}) { show_help(); exit 0; }
$opts{'m'} //= 'GET';
die "Only GET method is supported for database downloads.\n" if uc($opts{'m'}) ne 'GET';

make_path($opts{auth_dir}) unless -d $opts{auth_dir};
my $target_url = _resolve_target_url();
my $db_root    = _derive_db_root($target_url);

die "This OAuth client is for PubMLST (rest.pubmlst.org), but URL is: $target_url\n"
    unless $target_url =~ m{^https?://rest\.pubmlst\.org/}i;

main();

sub main {
    my ($session_token, $session_secret);
    my $data;
    my $retry_count = 0;

    while ($retry_count < 2) {
        ($session_token, $session_secret) = _retrieve_token('session_token');
        if (!defined $session_token || !defined $session_secret) {
            my $session_response = _get_session_token();
            ($session_token, $session_secret) = ($session_response->token, $session_response->token_secret);
        }

        my $response = _download_data($target_url, $session_token, $session_secret);
        if ($response->is_success) {
            $data = $response->content;
            last;
        } elsif ($response->code == 401) {
            warn "Session token expired or invalid. Renewing the global session token...\n";
            my $session_file = _token_path('session_token');
            unlink $session_file if -e $session_file;
            ($session_token, $session_secret) = (undef, undef);
            $retry_count++;
        } else {
            die "Download failed: " . $response->status_line . "\nURL: $target_url\n";
        }
    }

    die "Failed to download data after retry. Check PubMLST authorization/network.\n" unless defined $data;

    if ($opts{'o'}) {
        my $parent = dirname($opts{'o'});
        make_path($parent) if $parent && $parent ne '.' && !-d $parent;
        open(my $fh, '>', $opts{'o'}) or die "Cannot open output file $opts{'o'}: $!\n";
        binmode $fh;
        print $fh $data;
        close $fh;
    } else {
        print $data;
    }
}

sub _resolve_target_url {
    if (defined $opts{'u'} && $opts{'u'} ne '') {
        return $opts{'u'};
    }
    my $base = $opts{'b'} // '';
    my $route = $opts{'r'} // '';
    die "Provide --url URL, or both --base-url DB_ROOT and --route ROUTE.\n"
        unless $base ne '' && $route ne '';
    $base =~ s{/+$}{};
    $route =~ s{^/+}{};
    return "$base/$route";
}

sub _derive_db_root {
    my ($url) = @_;
    if ($url =~ m{^(https?://rest\.pubmlst\.org/db/[^/?#]+)}i) {
        return $1;
    }
    die "Cannot derive PubMLST database root from URL: $url\n";
}

sub _download_data {
    my ($url, $session_token, $session_secret) = @_;
    my %params = _parse_extra_params();

    my $request = Net::OAuth->request('protected resource')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $session_token,
        token_secret     => $session_secret,
        request_url      => $url,
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
        extra_params     => \%params,
    );
    $request->sign;
    die "OAuth signature verification failed!\n" unless $request->verify;

    my $ua = LWP::UserAgent->new(timeout => 120, agent => 'MLST-Typer-PubMLST-OAuth/2.0');
    return $ua->get($request->to_url);
}

sub _parse_extra_params {
    my %params;
    if ($opts{'a'}) {
        for my $pair (split /&/, $opts{'a'}) {
            my ($key, $value) = split /=/, $pair, 2;
            $params{$key} = defined $value ? $value : '';
        }
    }
    return %params;
}

sub _get_request_token {
    my $request = Net::OAuth->request('request token')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        request_url      => $db_root . '/oauth/get_request_token',
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
        callback         => 'oob'
    );
    $request->sign;
    die "COULDN'T VERIFY request-token signature.\n" unless $request->verify;

    say 'Getting request token...';
    my $ua = LWP::UserAgent->new(timeout => 120);
    my $res = $ua->request(GET $request->to_url);
    die "Failed to get request token: " . $res->status_line . "\n" unless $res->is_success;

    my $decoded_json = decode_json($res->content);
    my $request_response = Net::OAuth->response('request token')->from_hash($decoded_json);
    _write_token('request_token', $request_response->token, $request_response->token_secret);
    return $request_response;
}

sub _get_access_token {
    my ($request_token, $request_secret) = @_;
    unless ($request_token && $request_secret) {
        ($request_token, $request_secret) = _retrieve_token('request_token');
        unless ($request_token && $request_secret) {
            my $request_response = _get_request_token();
            ($request_token, $request_secret) = ($request_response->token, $request_response->token_secret);
        }
    }

    my $db_name = $db_root;
    $db_name =~ s{^.*/db/}{};
    my $authorize_url = "https://pubmlst.org/bigsdb?db=$db_name&page=authorizeClient&oauth_token=$request_token";
    say "\nPlease visit the following URL to authorize this client:";
    say $authorize_url;
    print "\nEnter verification code: ";
    my $verifier = <STDIN>;
    defined $verifier or die "No verification code received.\n";
    chomp $verifier;

    my $request = Net::OAuth->request('access token')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $request_token,
        token_secret     => $request_secret,
        verifier         => $verifier,
        request_url      => $db_root . '/oauth/get_access_token',
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
    );
    $request->sign;
    die "COULDN'T VERIFY access-token signature.\n" unless $request->verify;

    say "\nExchanging request token for global access token...";
    my $request_file = _token_path('request_token');
    unlink $request_file if -e $request_file;

    my $ua = LWP::UserAgent->new(timeout => 120);
    my $res = $ua->request(GET $request->to_url);
    die "Failed to get access token: " . $res->status_line . "\n" unless $res->is_success;

    my $decoded_json = decode_json($res->content);
    my $access_response = Net::OAuth->response('access token')->from_hash($decoded_json);
    _write_token('access_token', $access_response->token, $access_response->token_secret);
    say "Global access token saved in auth/access_token.";
    return $access_response;
}

sub _get_session_token {
    my ($access_token, $access_secret) = @_;
    unless ($access_token && $access_secret) {
        ($access_token, $access_secret) = _retrieve_token('access_token');
        unless ($access_token && $access_secret) {
            my $access_response = _get_access_token();
            ($access_token, $access_secret) = ($access_response->token, $access_response->token_secret);
        }
    }

    say "Requesting global session token...";
    my $request = Net::OAuth->request('protected resource')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $access_token,
        token_secret     => $access_secret,
        request_url      => $db_root . '/oauth/get_session_token',
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
    );
    $request->sign;
    die "COULDN'T VERIFY session-token signature.\n" unless $request->verify;

    my $ua = LWP::UserAgent->new(timeout => 120);
    my $res = $ua->request(GET $request->to_url);
    die "Failed to get session token: " . $res->status_line . "\n" unless $res->is_success;

    my $decoded_json = decode_json($res->content);
    my $session_response = Net::OAuth->response('access token')->from_hash($decoded_json);
    _write_token('session_token', $session_response->token, $session_response->token_secret);
    say "Global session token saved in auth/session_token.";
    return $session_response;
}

sub _token_path {
    my ($token_name) = @_;
    return "$opts{auth_dir}/$token_name";
}

sub _retrieve_token {
    my ($token_name) = @_;
    my $path = _token_path($token_name);
    return unless -e $path;
    my $config = Config::Tiny->read($path) or die "Failed to read token file: $path\n";
    return ($config->{_}->{token}, $config->{_}->{secret});
}

sub _write_token {
    my ($token_name, $token, $secret) = @_;
    my $path = _token_path($token_name);
    my $config = Config::Tiny->new;
    $config->{_} = { token => $token, secret => $secret };
    $config->write($path) or die "Failed to write token file: $path\n";
    chmod 0600, $path;
    return 1;
}

sub show_help {
    print <<'HELP';
PubMLST authenticated download tool (global OAuth tokens)

Usage:
  perl RESTful_API_download_pubmlst.pl --url FULL_URL --output FILE
  perl RESTful_API_download_pubmlst.pl --base-url DB_ROOT --route ROUTE --output FILE

Examples:
  perl RESTful_API_download_pubmlst.pl \
    --url https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/aroE/alleles_fasta \
    --output aroE.fas

Authentication behavior:
  - ONE auth/access_token is reused for all PubMLST databases/schemes.
  - ONE auth/session_token is reused for all PubMLST databases/schemes.
  - If session_token is missing/expired, it is renewed from access_token.
  - If access_token is missing, the original request-token -> browser authorization
    -> access-token workflow is started once.

Options:
  --url, -u URL        Complete PubMLST REST resource URL.
  --base-url, -b URL   Database root, e.g. https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef
  --route, -r ROUTE    Relative API route when --base-url is used.
  --output, -o FILE    Save response to file.
  --arguments, -a STR  Optional query parameters.
  --auth-dir DIR       Global token directory (default: ./auth next to this script).
HELP
}
